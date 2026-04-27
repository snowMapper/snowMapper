"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description:
Builds the classifier used for the reconstruction of snow cover, using in-situ 
measurements of snow depth that have been previously preprocessed into binary snow cover.

Input parameters:
- training_ee (ee.FeatureCollection): Training dataset of in-situ measurements (points),
                                      that includes binary snow cover, as well as all
                                      metadata bands included in the image collection
                                      for reconstruction to be performed.
- settings (dict): Settings for the selected classifier.
- active_classifier (str): Classifier that has been selected in the configuration file 
                           for the reconstruction.
- classifier_path (str): Path where EE classifier asset will be stored.

Internal functions:
- register_method(): Finds the classifier that needs to be applied.
- train_classifier(): Applies the selected classifier.

Output:
- ee.Classifier used in .reconstruct_daily() / .reconstruct_monthly().
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import ee
import os
import pandas as pd
import matplotlib.pyplot as plt

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Method registry
#===============================================================================
CLASSIFIER_REGISTRY = {}

def register_classifier(name):
    def decorator(func):
        CLASSIFIER_REGISTRY[name] = func
        return func
    return decorator

#===============================================================================
# Method: Random Forest
#===============================================================================
@register_classifier("random_forest")
def classifier_random_forest(training_ee, settings, input_vars):
    trained_classifier = (
        ee.Classifier
          .smileRandomForest(
              settings['numberOfTrees'],
              settings['variablesPerSplit'],
              settings['minLeafPopulation'],
              settings['bagFraction'],
              settings['maxNodes'],
              settings['seed']
          )
          .setOutputMode(settings['outputMode'])
          .train(
              features = training_ee,
              classProperty = settings['classProperty'],
              inputProperties = input_vars
          )
    )
    return trained_classifier

#===============================================================================
# Method: Gradient Tree Boost
#===============================================================================
@register_classifier("gradient_tree_boost")
def classifier_gradient_tree_boost(training_ee, settings, input_vars):
    trained_classifier = (
        ee.Classifier
          .smileGradientTreeBoost(
              settings['numberOfTrees'],
              settings['shrinkage'],
              settings['samplingRate'],
              settings['maxNodes'],
              settings['loss'],
              settings['seed']
          )
          .setOutputMode(settings['outputMode'])
          .train(
              features = training_ee,
              classProperty = settings['classProperty'],
              inputProperties = input_vars
          )
    )
    return trained_classifier

#===============================================================================
# Method: Minimum Distance
#===============================================================================
@register_classifier("minimum_distance")
def classifier_minimum_distance(training_ee, settings, input_vars):
    trained_classifier = (
        ee.Classifier
          .minimumDistance(
              settings['metric'], 
              settings['kNearest']
          )
          .setOutputMode(settings['outputMode'])
          .train(
              features = training_ee,
              classProperty = settings['classProperty'],
              inputProperties = input_vars
          )
    )
    return trained_classifier

#===============================================================================
# Method: Classification and Regression Trees (CART)
#===============================================================================
@register_classifier("cart")
def classifier_cart(training_ee, settings, input_vars):
    trained_classifier = (
        ee.Classifier \
          .smileCart(
              settings['maxNodes'], 
              settings['minLeafPopulation']
          )
          .setOutputMode(settings['outputMode'])
          .train(
              features = training_ee,
              classProperty = settings['classProperty'],
              inputProperties = input_vars
          )
    )
    return trained_classifier

#===============================================================================
# Method: k-Nearest Neighbor (k-NN) 
#===============================================================================
@register_classifier("k_nearest_neighbor")
def classifier_k_nearest_neighbor(training_ee, settings, input_vars):
    trained_classifier = (
        ee.Classifier
          .smileKNN(
              settings['k'],
              settings['searchMethod'],
              settings['metric']
          )
          .setOutputMode(settings['outputMode'])
          .train(
              features = training_ee,
              classProperty = settings['classProperty'],
              inputProperties = input_vars
          )
    )
    return trained_classifier

#===============================================================================
# Method: Naive Bayes
#===============================================================================
@register_classifier("naive_bayes")
def classifier_naive_bayes(training_ee, settings, input_vars):
    trained_classifier = (
        ee.Classifier
          .smileNaiveBayes(
              settings['lambda']
          )
          .setOutputMode(settings['outputMode'])
          .train(
              features = training_ee,
              classProperty = settings['classProperty'],
              inputProperties = input_vars
          )
    )
    return trained_classifier

#===============================================================================
# Method: Support Vector Machine
#===============================================================================
@register_classifier("support_vector_machine")
def classifier_support_vector_machine(training_ee, settings, input_vars):
    trained_classifier = (
        ee.Classifier
          .libsvm(
              settings['decisionProcedure'],
              settings['svmType'],
              settings['kernelType'],
              settings['shrinking'],
              settings['degree'],
              settings['gamma'],
              settings['coef0'],
              settings['cost'],
              settings['nu'],
              settings['terminationEpsilon'],
              settings['lossEpsilon'], 
              settings['oneClass']
          )
          .setOutputMode(settings['outputMode'])
          .train(
              features = training_ee,
              classProperty = settings['classProperty'],
              inputProperties = input_vars
          )
    )
    return trained_classifier

#===============================================================================
# Train classifier
#===============================================================================
def train_classifier(training_ee, settings, input_vars, active_classifier, ee_project, output_root):

    classifier = CLASSIFIER_REGISTRY[active_classifier](training_ee, settings, input_vars)

    asset_id = f"projects/{ee_project}/assets/snowMapper/{active_classifier}_classifier"
    task = ee.batch.Export.classifier.toAsset(
        classifier=classifier,
        description='classifier_export',
        assetId=asset_id
    )
    task.start()

    # Feature Importance
    importance_dict = classifier.explain().get('importance').getInfo()
    total = sum(importance_dict.values())
    normalized_dict = {k: (v / total) * 100 for k, v in importance_dict.items()}
    sorted_items = sorted(normalized_dict.items(), key=lambda x: x[1], reverse=True)
    features, importances = zip(*sorted_items)
    
    df = pd.DataFrame({'features':features, 'importance':importances})
    
    # Plot feature importance
    plt.figure(figsize=(10, 6))
    plt.barh(df["features"], df["importance"])
    plt.xlabel('Normalised Importance (%)')
    
    plt.gca().invert_yaxis()  # Highest importance at top
    plt.tight_layout()  
    fig_path = os.path.join(output_root, "ml_feature_importance.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    
    return classifier