#include "../rtmonopitch.h"

#include "util_p.hpp"
#include "../rtpyin.h"
#include "../rthmm.h"
#include <vector>

using namespace ReVoice;

typedef struct RvRTMonoPitchProcessor
{
  RvRTMonoPitchProcessorParameter param;
  RvRTSparseHMM *hmmModel;
  RvReal *obsTemp;
  int *decodeTemp;
  std::vector<std::vector<std::pair<RvReal, RvReal>>> obsProbList;
  std::vector<bool> obsSilentList;
  int nState, nTrans;

} RvRTMonoPitchProcessor;

RvRTMonoPitchProcessorParameter *rvCreateRTMonoPitchProcessorParameter(int hopSize, RvReal samprate, int nSemitone, RvReal maxTransSemitone, RvReal minFreq)
{
  RvReal nyq = samprate / 2.0;
  rvAssert(hopSize > 0, "invalid hopSize");
  rvAssert(samprate > 0.0, "invalid samprate");
  rvAssert(nSemitone > 0, "invalid nSemitone");
  rvAssert(maxTransSemitone > 0.0, "invalid maxTransSemitone");
  rvAssert(minFreq > 0.0 && minFreq < nyq, "invalid minFreq");

  auto p = new RvRTMonoPitchProcessorParameter;
  p->hopSize = hopSize;
  p->samprate = samprate;
  p->nSemitone = nSemitone;
  p->maxTransSemitone = maxTransSemitone;
  p->minFreq = minFreq;

  p->binPerSemitone = 5;
  p->transSelf = 0.999;
  p->yinTrust = 0.5;
  p->energyThreshold = 1e-8;
  p->maxObsLength = 128;

  return p;
}

RvRTMonoPitchProcessorParameter *rvCreateRTMonoPitchProcessorParameterFromRTPYin(const RvRTPYinProcessorParameter *param)
{
  int hopSize = param->hopSize;
  RvReal samprate = param->samprate;
  int nSemitone = static_cast<int>(std::ceil(std::log2(param->maxFreq / param->minFreq) * 12.0));
  RvReal maxTransSemitone = (param->hopSize / param->samprate) / (256.0 / 44100.0) * 3.0;
  RvReal minFreq = param->minFreq;
  return rvCreateRTMonoPitchProcessorParameter(hopSize, samprate, nSemitone, maxTransSemitone, minFreq);
}

void rvDestroyRTMonoPitchProcessorParameter(RvRTMonoPitchProcessorParameter *param)
{
  rvAssert(param, "param cannot be nullptr");
  delete param;
}

RvRTMonoPitchProcessor *rvCreateRTMonoPitchProcessor(const RvRTMonoPitchProcessorParameter *param)
{
  rvAssert(param, "param cannot be nullptr");
  RvReal nyq = param->samprate / 2.0;
  rvAssert(param->hopSize > 0, "invalid hopSize");
  rvAssert(param->samprate > 0.0, "invalid samprate");
  rvAssert(param->nSemitone > 0, "invalid nSemitone");
  rvAssert(param->maxTransSemitone > 0.0, "invalid maxTransSemitone");
  rvAssert(param->minFreq > 0.0 && param->minFreq < nyq, "invalid minFreq");
  rvAssert(param->binPerSemitone > 0, "invalid binPerSemitone");
  rvAssert(param->transSelf >= 0.0 && param->transSelf <= 1.0, "invalid transSelf");
  rvAssert(param->yinTrust >= 0.0 && param->yinTrust <= 1.0, "invalid yinTrust");
  rvAssert(param->energyThreshold >= 0.0, "invalid energyThreshold");
  rvAssert(param->maxObsLength > 0, "invalid maxObsLength");

  auto self = new RvRTMonoPitchProcessor;
  self->param = *param;

  // create model
  {
    int nBin = param->nSemitone * param->binPerSemitone;
    int halfMaxTransBin = static_cast<int>(std::round(param->maxTransSemitone * static_cast<RvReal>(param->binPerSemitone) / 2.0));
    int nState = 2 * nBin;
    int nTrans = 4 * (nBin * (2 * halfMaxTransBin + 1) - halfMaxTransBin * (halfMaxTransBin + 1));

    auto init = RVALLOC(RvReal, nState);
    auto frm = RVALLOC(int, nTrans);
    auto to = RVALLOC(int, nTrans);
    auto transProb = RVALLOC(RvReal, nTrans);

    std::fill(init, init + nState, 1.0 / static_cast<RvReal>(nState));
    int iA = 0;
    for(int iBin = 0; iBin < nBin; ++iBin)
    {
      int theoreticalMinNextBin = iBin - halfMaxTransBin;
      int minNextBin = std::max(iBin - halfMaxTransBin, 0);
      int maxNextBin = std::min(iBin + halfMaxTransBin, nBin - 1);
      std::vector<RvReal> weights(maxNextBin - minNextBin + 1, 0.0);
      for(int i = minNextBin; i < maxNextBin + 1; ++i)
        weights[i - minNextBin] = i <= iBin ? (i - theoreticalMinNextBin + 1.0) : (iBin - theoreticalMinNextBin + 1.0 - (i - iBin));
      RvReal weightSum = std::accumulate(weights.cbegin(), weights.cend(), 0.0);

      // trans to close pitch
      for(int i = minNextBin; i < maxNextBin + 1; ++i)
      {
        frm[iA] = iBin;
        to[iA] = i;
        transProb[iA] = weights[i - minNextBin] / weightSum * param->transSelf;

        frm[iA + 1] = iBin;
        to[iA + 1] = i + nBin;
        transProb[iA + 1] = weights[i - minNextBin] / weightSum * (1.0 - param->transSelf);

        frm[iA + 2] = iBin + nBin;
        to[iA + 2] = i + nBin;
        transProb[iA + 2] = weights[i - minNextBin] / weightSum * param->transSelf;

        frm[iA + 3] = iBin + nBin;
        to[iA + 3] = i;
        transProb[iA + 3] = weights[i - minNextBin] / weightSum * (1.0 - param->transSelf);
        iA += 4;
      }
    }

    self->hmmModel = rvCreateRTSparseHMM(init, frm, to, transProb, nState, nTrans, param->maxObsLength);
    self->nState = nState;
    self->nTrans = nTrans;

    rvFree(transProb);
    rvFree(to);
    rvFree(frm);
    rvFree(init);
  }

  self->obsTemp = RVALLOC(RvReal, self->nState);
  self->decodeTemp = RVALLOC(int, param->maxObsLength);
  return self;
}

const RvRTMonoPitchProcessorParameter *rvRTMonoPitchParam(const RvRTMonoPitchProcessor *self)
{ return &(self->param); }

int rvMonoPitchNextOutputLength(const RvRTMonoPitchProcessor *self)
{ return std::min(static_cast<int>(self->obsProbList.size() + 1), self->param.maxObsLength); }

void rvMonoPitchDumpObsTemp(const RvRTMonoPitchProcessor *self, RvReal *out)
{ std::copy(self->obsTemp, self->obsTemp + self->nState, out); }

int rvCallRTMonoPitch(RvRTMonoPitchProcessor *self, const RvReal *x, const RvReal *obsProb, int nObsProb, RvReal *out)
{
  rvAssert(x, "x cannot be nullptr");
  rvAssert(obsProb || nObsProb == 0, "obsProb cannot be nullptr with non-zero nObsProb");
  rvAssert(nObsProb >= 0, "invalid nObsProb");
  rvAssert(out, "out cannot be nullptr");

  auto &p = self->param;
  int nBin = p.nSemitone * p.binPerSemitone;
  int nState = self->nState;
  RvReal maxFreq = p.minFreq * std::pow(2.0, static_cast<RvReal>(p.nSemitone) / 12.0);
  // calc state prob
  {
    std::fill(self->obsTemp, self->obsTemp + self->nState, 0.0);
    RvReal probYinPitched = 0.0;
    for(int i = 0; i < nObsProb; ++i)
    {
      RvReal freq = obsProb[i * 2];
      RvReal prob = obsProb[i * 2 + 1];
      if(freq < p.minFreq || freq > maxFreq)
      {
        if(freq <= 0.0)
          break;
        continue;
      }
      int iBin = clip(0, static_cast<int>(std::round(std::log2(freq / p.minFreq) * 12.0 * static_cast<RvReal>(p.binPerSemitone))), nBin - 1);
      self->obsTemp[iBin] = prob;
      probYinPitched += prob;
    }
    RvReal probReallyPitched = p.yinTrust * probYinPitched;
    if(probYinPitched > 0.0)
    {
      RvReal v = probReallyPitched / probYinPitched;
      for(int i = 0; i < nBin; ++i)
        self->obsTemp[i] *= v;
    }
    RvReal v = (1.0 - probReallyPitched) / static_cast<RvReal>(nBin);
    for(int i = nBin; i < nState; ++i)
      self->obsTemp[i] = v;
    for(int i = 0; i < nState; ++i)
      self->obsTemp[i] = std::max(0.0, self->obsTemp[i]) + 1e-5;
  }

  // feed and decode
  int currObsLength = std::min(static_cast<int>(self->obsProbList.size() + 1), p.maxObsLength);
  rvRTSparseHMMFeed(self->hmmModel, self->obsTemp);
  
  int nDecoded = rvRTSparseHMMViterbiDecode(self->hmmModel, self->decodeTemp, currObsLength);
  rvAssert(nDecoded == currObsLength, "internal error");

  // save state
  bool isSilent;
  {
    RvReal dcOffset = mean(x, p.hopSize * 2);
    RvReal meanEnergy = 0.0;
    for(int i = 0; i < p.hopSize * 2; ++i)
    {
      RvReal v = x[i] - dcOffset;
      meanEnergy += v * v;
    }
    meanEnergy /= static_cast<RvReal>(p.hopSize * 2);
    isSilent = meanEnergy < p.energyThreshold;
  }
  if(self->obsProbList.size() >= p.maxObsLength)
  {
    self->obsProbList.erase(self->obsProbList.begin());
    self->obsSilentList.erase(self->obsSilentList.begin());
  }
  rvAssert(self->obsProbList.size() < p.maxObsLength && self->obsProbList.size() == self->obsSilentList.size(), "internal error");

  {
    std::vector<std::pair<RvReal, RvReal>> obsProbVec(nObsProb);
    for(int i = 0; i < nObsProb; ++i)
      obsProbVec[i] = std::make_pair(obsProb[i * 2], obsProb[i * 2 + 1]);
    self->obsProbList.push_back(std::move(obsProbVec));
    self->obsSilentList.push_back(isSilent);
  }

  // extract frequency from path
  for(int iHop = 0; iHop < currObsLength; ++iHop)
  {
    RvReal bestFreq;
    if(self->decodeTemp[iHop] < nBin)
    {
      RvReal hmmFreq = p.minFreq * std::pow(2.0, static_cast<RvReal>(self->decodeTemp[iHop]) / (12.0 * static_cast<RvReal>(p.binPerSemitone)));
      if(self->obsProbList[iHop].size() == 0)
        bestFreq = hmmFreq;
      else
      {
        int iNearest = -1;
        RvReal nearestDistance = std::numeric_limits<RvReal>::infinity();
        int n = static_cast<int>(self->obsProbList[iHop].size());
        for(int i = 0; i < n; ++i)
        {
          RvReal distance = std::abs(self->obsProbList[iHop][i].first - hmmFreq);
          if(distance < nearestDistance)
          {
            nearestDistance = distance;
            iNearest = i;
          }
        }
        rvAssert(iNearest >= 0, "internal error");
        bestFreq = self->obsProbList[iHop][iNearest].first;
        if(bestFreq < p.minFreq || bestFreq > maxFreq || std::abs(std::log2(bestFreq / p.minFreq) * 12.0 * static_cast<RvReal>(p.binPerSemitone) - self->decodeTemp[iHop] > 1.0))
          bestFreq = hmmFreq;
      }
    }
    else
      bestFreq = -p.minFreq * std::pow(2.0, (static_cast<RvReal>(self->decodeTemp[iHop] - nBin)) / (12.0 * static_cast<RvReal>(p.binPerSemitone)));
    out[iHop] = bestFreq;
  }

  // mark unvoiced->voiced bound as voiced
  for(int iHop = 1; iHop < currObsLength; ++iHop)
  {
    if(out[iHop - 1] <= 0.0 && out[iHop] > 0.0)
    {
      int windowSize = std::max(static_cast<int>(std::ceil(p.samprate / out[iHop] * 4.0)), p.hopSize * 2);
      if(windowSize % 2 != 0)
        windowSize += 1;
      int frameOffset = static_cast<int>(std::round(static_cast<RvReal>(windowSize) / static_cast<RvReal>(p.hopSize * 2)));
      for(int i = std::max(0, iHop - frameOffset); i < iHop; ++i)
        out[i] = out[iHop];
    }
  }

  // mark silent frame as unvoiced
  for(int iHop = 0; iHop < currObsLength; ++iHop)
  {
    if(out[iHop] > 0.0 && self->obsSilentList[iHop])
      out[iHop] = 0.0;
  }

  return currObsLength;
}

void rvDestroyRTMonoPitchProcessor(RvRTMonoPitchProcessor *self)
{
  rvFree(self->obsTemp);
  rvFree(self->decodeTemp);
  rvDestroyRTSparseHMM(self->hmmModel);
}