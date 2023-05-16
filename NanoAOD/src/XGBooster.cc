#include "Bmm5/NanoAOD/interface/XGBooster.h"
#include <assert.h>
#include <math.h>
#include <stdexcept>
#include <algorithm>

XGBooster::XGBooster(std::string model_file)
{
  int status = XGBoosterCreate(NULL, 0, &booster_);
  assert(status==0 && "Failed to create XGBooster");
  status = XGBoosterLoadModel(booster_, model_file.c_str());
  assert(status==0 && "Failed to load XGBoost model");
  XGBoosterSetParam(booster_, "nthread", "1");
};

void XGBooster::reset()
{
  std::fill(features_.begin(), features_.end(), std::nan(""));
}


void XGBooster::addFeature(std::string name){
  features_.push_back(0);
  feature_name_to_index_[name] = features_.size()-1;
}

void XGBooster::set(std::string name, float value){
  features_.at(feature_name_to_index_[name]) = value;
}

float XGBooster::predict()
{
  float result(-999.);

  // check if all feature values are set properly
  for (auto feature: features_)
    if (isnan(feature)){
      reset();
      return result;
    }
  
  DMatrixHandle dvalues;
  XGDMatrixCreateFromMat(&features_[0], 1, features_.size(), 0., &dvalues);
    
  bst_ulong out_len=0;
  const float* score;

  auto ret = XGBoosterPredict(booster_, dvalues, 0, 0, 0, &out_len, &score);

  XGDMatrixFree(dvalues);

  if (ret==0) {
    assert(out_len==1 && "Unexpected prediction format");
    result = score[0];
  }
  
  reset();
  
  return result;    
}
