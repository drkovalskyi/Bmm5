#include <memory>
#include <string>
#include <vector>
#include <map>
#include <xgboost/c_api.h>

class XGBooster
{       
public:
    XGBooster(std::string model_file);

    /// Features need to be entered in the order they are used
    /// in the model
    void addFeature(std::string name);
    
    void set(std::string name, float value);
    
    float predict();

private:
    std::vector<float> features_;
    std::map<std::string,unsigned int> feature_name_to_index_;
    BoosterHandle  booster_;
};
