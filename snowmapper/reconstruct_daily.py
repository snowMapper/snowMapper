"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
Full daily reconstruction of snow cover, using a machine learning classifier trained on 
binarised snow depth in situ observations (snow cover). It is meant to run sequentially as
part of an iteration - i.e. using iterate() instead of map() - since it requires the 
previous day's snow cover state to inform the current day's variables and reconstruction.

Input parameters:
- collection (ee.ImageCollection): Images with metadata.
- initialisation_date (str): Initialisation ('%Y-%m-%d') 1 day before the start-date.
- start_date (str): Start-date of the season ('%Y-%m-%d').
- classifier (classifier): The random forest classifier developed in sm.train_classifier().
- decision_tree_settings (dict): Contains booleans and thresholds for decision tree-based 
                                 gap-filling.
- input_vars (list): The input variables to be used for the reconstruction.
- output_vars (list): The output variables to be saved after the reconstruction.

Internal functions:
- get_prev_day_vars(): Add the previous day's cumulative degree-day variables
- dt_gapfill(): Decision tree-based gap-filing using the previous day's snow cover,
                current day's temperature & precipitation, current months multi-year
                monthly snow cover probability, and MODIS snow cover data.
- ml_gapfill(): Machine learning-based spatial gap-filling using a classifier.
- evaluate(): Produce evaluation bands true positive (tp), true negative (tn), false
              positive (fp), and false negative (fn). These are produced only for pixels
              where observations exist (not for sc_obs==2).
- assimilate(): Assimilate observations after img reconstruction.
- update_sum_vars(): Update cumulative degree-day variables post-reconstruction.

Output:
- ee.ImageCollection of reconstructed images.
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import ee

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Reconstruct all binary snow images
#===============================================================================
def reconstruct_daily(collection, initialisation_date, start_date, decision_tree_settings, classifier, input_vars, output_vars):

    initialisation_img = ee.Image(collection.first())
    post_initialisation_col = ee.ImageCollection(collection.filter(ee.Filter.date(initialisation_date, start_date).Not()))
    # post_initialisation_col = ee.ImageCollection(collection.filter(ee.Filter.calendarRange(START_MONTH, START_MONTH, 'month').Not()))

    #------------------------------------------
    # Prepare for reconstruction iteration
    #------------------------------------------
    def reconstruct_img(img, prev_imgs, decision_tree_settings, classifier, input_vars):

        img = ee.Image(img)
        acquisition_date = img.get('system:time_start')
        prev_imgs = ee.List(prev_imgs)
        prev_img = ee.Image(prev_imgs.get(-1))

        #------------------------------------------
        # Degree-day parameters
        #------------------------------------------
        def get_prev_day_vars(img, prev_img):
            sc_sum_prev = prev_img.select(['sc_sum'], ['sc_sum_prev'])
            hd_sum_prev = prev_img.select(['hd_sum'], ['hd_sum_prev'])
            cd_sum_prev = prev_img.select(['cd_sum'], ['cd_sum_prev'])
            hd_precip_sum_prev = prev_img.select(['hd_precip_sum'], ['hd_precip_sum_prev'])
            cd_precip_sum_prev = prev_img.select(['cd_precip_sum'], ['cd_precip_sum_prev'])

            return ee.Image(
                img.addBands([
                    sc_sum_prev, hd_sum_prev, cd_sum_prev, hd_precip_sum_prev, cd_precip_sum_prev
                ], overwrite=True)
            )

        #------------------------------------------
        # Decision tree-based gapfilling
        #------------------------------------------
        def dt_gapfill(img, decision_tree_settings):
            # Mask no-data pixels (sc == 2)
            gap_mask = img.select('sc').eq(2)

            # Decision tree 1: Rule based on sc_sum_prev, t2mHres, precipHres
            if decision_tree_settings["temperature"]:
                cond1_mask_1 = gap_mask.And(
                    img.select('sc_sum_prev').gte(1)
                       .And(img.select('t2mHres').lte(decision_tree_settings["FREEZING_TEMP"]))
                )
                cond1_mask_0 = gap_mask.And(
                    img.select('sc_sum_prev').eq(0)
                       .And(img.select('t2mHres').gt(decision_tree_settings["FREEZING_TEMP"]))
                       .And(img.select('precipHres').eq(0))
                )
    
                # Apply decision tree 1
                sc_filled = (
                    img.select('sc')
                       .where(cond1_mask_1, 1)
                       .where(cond1_mask_0, 0)
                )
    
                # Update gap_mask
                gap_mask = sc_filled.eq(2)

            # Decision tree 2: Use sc_modis to fill gaps
            if decision_tree_settings["modis"]:
                cond2_mask_1 = gap_mask.And(img.select('sc_modis').eq(1))
                cond2_mask_0 = gap_mask.And(img.select('sc_modis').eq(0))
    
                # Apply decision tree 2
                sc_filled = (
                    sc_filled.where(cond2_mask_1, 1)
                             .where(cond2_mask_0, 0)
                )
    
                # Update gap_mask
                gap_mask = sc_filled.eq(2)

            # Decision tree 3: Use sc_probab to fill based on snow probability threshold
            if decision_tree_settings["sc_probab"]:
                cond3_mask_1 = gap_mask.And(img.select('sc_probab').gte(decision_tree_settings["SNOW_PROBAB_THRES"]))
                cond3_mask_0 = gap_mask.And(img.select('sc_probab').lte(1 - decision_tree_settings["SNOW_PROBAB_THRES"]))
    
                # Apply decision tree 3
                sc_filled = (
                    sc_filled.where(cond3_mask_1, 1)
                             .where(cond3_mask_0, 0)
                )
            
            # Define mask of where we successfully filled gaps
            filled_mask = img.select('flag_dt').eq(1).And(sc_filled.neq(2))
            
            # Update flag_dt: set to 0 where previously 1 and now filled
            flag_dt_updated = img.select('flag_dt').where(filled_mask, 0)
            
            # Add bands: update 'sc' and 'flag_dt'
            img_filled = img.addBands(
                [sc_filled.rename('sc').toUint8(),
                 flag_dt_updated.rename('flag_dt').toUint8()
                ],
                overwrite=True
            )
            return img_filled
            
        #------------------------------------------
        # Machine learning-based gapfill
        #------------------------------------------
        def ml_gapfill(img, classifier, input_vars):
            # Classify only where sc == 2
            gap_mask = img.select('sc').eq(2)
            gap_img = img.updateMask(gap_mask).select(input_vars)
            classified_gap = gap_img.classify(classifier)

            # Replace sc == 2 pixels with classified values
            original_sc = img.select('sc')
            filled_sc = original_sc.where(gap_mask, classified_gap)

            # Return updated image (replacing the 'sc' band)
            output = img.addBands(filled_sc.rename('sc'), overwrite=True)
            return output

        #------------------------------------------
        # Evaluate the reconstruction performance
        #------------------------------------------
        def evaluate(img):
            sc = img.select('sc')
            sc_obs = img.select('sc_obs')

            # Create a mask where sc_obs is not 2 (i.e., valid data)
            valid_mask = sc_obs.neq(2)

            # Predicted snow and observed snow
            predicted = sc.updateMask(valid_mask)
            observed = sc_obs.updateMask(valid_mask)

            # True Positive: predicted = 1 and observed = 1
            tp = predicted.eq(1).And(observed.eq(1)).toUint8().rename('tp')

            # True Negative: predicted = 0 and observed = 0
            tn = predicted.eq(0).And(observed.eq(0)).toUint8().rename('tn')

            # False Positive: predicted = 1 and observed = 0
            fp = predicted.eq(1).And(observed.eq(0)).toUint8().rename('fp')

            # False Negative: predicted = 0 and observed = 1
            fn = predicted.eq(0).And(observed.eq(1)).toUint8().rename('fn')

            # Return original image with new bands
            return img.addBands([tp, tn, fp, fn], overwrite=True)

        #------------------------------------------
        # Replace 'sc' band with available 'sc_obs' band
        #------------------------------------------
        def assimilate(img):
            sc = img.select('sc')
            sc_obs = img.select('sc_obs')

            # Create a mask where sc_obs is not equal to 2
            valid_mask = sc_obs.neq(2)

            # Replace sc values where valid_mask is true with sc_obs values
            filled_sc = sc.where(valid_mask, sc_obs)

            # Replace 'sc' band in the original image with the filled one
            assimilated_img = img.addBands(filled_sc.rename('sc'), overwrite=True)
            return assimilated_img

        #------------------------------------------
        # Update cumulative snow cover value
        #------------------------------------------
        def update_sum_vars(img):
            sc = img.select('sc')
            hd = img.select('hd')
            cd = img.select('cd')
            hd_precip = img.select('hd_precip')
            cd_precip = img.select('cd_precip')
            sc_sum_prev = img.select('sc_sum_prev')
            hd_sum_prev = img.select('hd_sum_prev')
            cd_sum_prev = img.select('cd_sum_prev')
            hd_precip_sum_prev = img.select('hd_precip_sum_prev')
            cd_precip_sum_prev = img.select('cd_precip_sum_prev')

            # Compute new bands
            sc_sum = (sc_sum_prev.add(sc)).multiply(sc).toUint16().rename('sc_sum')
            hd_sum = (hd_sum_prev.add(hd)).multiply(sc).toFloat().rename('hd_sum')
            cd_sum = (cd_sum_prev.add(cd)).multiply(sc).toFloat().rename('cd_sum')    
            hd_precip_sum = (hd_precip_sum_prev.add(hd_precip)).multiply(sc).toFloat().rename('hd_precip_sum')
            cd_precip_sum = (cd_precip_sum_prev.add(cd_precip)).multiply(sc).toFloat().rename('cd_precip_sum')    

            return ee.Image(
                img.addBands([sc_sum, hd_sum, cd_sum, hd_precip_sum, cd_precip_sum], overwrite=True)
            )

        #------------------------------------------
        # Apply the above
        #------------------------------------------
        updated_vars_img = ee.Image(get_prev_day_vars(img, prev_img))
        dt_gapfilled = ee.Image(dt_gapfill(updated_vars_img, decision_tree_settings))
        gapfilled_img = ee.Image(ml_gapfill(dt_gapfilled, classifier, input_vars))
        eval_img = ee.Image(evaluate(gapfilled_img))
        assimilated_img = ee.Image(assimilate(eval_img))
        final_img = ee.Image(update_sum_vars(assimilated_img).set('system:time_start', acquisition_date))
        return ee.List(prev_imgs).add(final_img)

    #------------------------------------------
    # Run iteration
    #------------------------------------------
    reconstructed_col = ee.List(
        post_initialisation_col.iterate(
            lambda img, initialisation_img:
            reconstruct_img(img, initialisation_img, decision_tree_settings, classifier, input_vars),  # This is the function to pass to iterate().
            ee.List([initialisation_img])                                                              # This is the initial state image.
        )
    )

    reconstructed_col = ee.ImageCollection(reconstructed_col).map(
        lambda img: img.select(output_vars)
    )

    return reconstructed_col.sort('system:time_start')