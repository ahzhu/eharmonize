#!/usr/bin/env python
import eharmonize.funcs.utils as efu 
import pandas as pd
import numpy as np
import click
import warnings
# from datetime import datetime
import sys
import os
import json

warnings.filterwarnings("ignore")

@click.group(context_settings={'help_option_names': ['-h', '--help']})
def eharmonize():
    '''
    eharmonize (ENIGMA harmonization) is a python-based tool for harmonizing outputs of the ENIGMA-DTI pipeline to a provided reference

    Input CSV files should include the following columns: subjectID, Age, Sex (coded 0 for females, 1 for males), SITE (optional; if multi-site), Dx (optional; coded "case" and "control"), ENIGMA-DTI outputs
    '''
    pass

@eharmonize.command()
@click.option("--incsv", metavar="CSV", help="File path of the CSV file with the necessary covariates and FA measures")
@click.option("--outdir", metavar="DIR", help="Directory to write out outputs to")
# @click.option("--site", default="SITE", metavar="column", help="(optional) Column name to indicate site/study/protocol to harmonize by"
@click.option("--reference", show_default=True, default='v0.0', metavar="version", help="(optional) Version number of desired reference")
@click.option("--rerun", is_flag=True, help="(optional) if used, will overwrite previous outputs")
def harmonize_FA(incsv, outdir, reference, rerun):
    
    # Settings Check

    if not os.path.exists(incsv):
        raise OSError("--incsv input does not exist:\n%s" %incsv)

    if os.path.exists(os.path.join(outdir, "harmonized_FA_%s.model" %reference)):
        if rerun:
            os.remove(os.path.join(outdir, "harmonized_FA_%s.model" %reference))
        else:
            raise ValueError("Model file already exists. Please use --rerun flag if you mean to overwrite")

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # now = datetime.now()
    # dt_string = now.strftime("%B %d, %Y %H:%M:%S")
    global txtlog
    txtlog = efu.logstr()
    log = txtlog.run_settings
    log["reference"] = reference

    dfI = pd.read_csv(incsv)
    covars = ["Age", "Sex", "SITE"]
    if "SITE" not in dfI.columns:
        dfI.loc[:, "SITE"] = "enigma"
        site_msg = "\nNo 'SITE' column provided, so we're going to assume that all data comes from the same site.\n"
        txtlog.stdO_file(site_msg)

    if "DX" in dfI.columns:
        dfI.rename(columns={"DX": "Dx"}, inplace=True)
    if "Dx" in dfI.columns:
        dfI.loc[dfI.Dx.notnull(), "Dx"] = dfI.loc[dfI.Dx.notnull(), "Dx"].str.strip().str.lower()
    cohort_info = efu.input_check(dfI, covars[:-1])
    log.update(cohort_info)

    log["FA"] = {}
    pkg_settings, dfR = efu.load_version(reference, "FA")
    roi_columns = dfR.columns[3:] 
    dfF, roi_dict, outM = efu.column_match(dfI, roi_columns, "FA")
    txtlog.stdO_file(outM)
    txtlog.stdO_file("\nUsing ROIS:\n%s\n" % ", ".join(roi_dict.keys()))
    log["FA"]["ROIs"] = list(roi_dict.keys())

    # Subject Check

    age_excludes, outM = efu.age_check(dfF, pkg_settings)
    txtlog.stdO_file(outM)
    dfF = dfF.loc[~age_excludes]

    # Preprocessing 

    if "Dx" in dfF.columns:
        controlDF = dfF.loc[dfF.Dx== "control"]
        if controlDF.empty:
            txtlog.stdO_file("\nWe appear to have detected a case only input.\nWe're not set up at the moment, so if that's the case, sorry.\nIf not, please check inputs and try again.\n")
            txtlog.to_file(os.path.join(outdir, "harmonized_FA_%s.txt" %reference))
            sys.exit(1)
        caseDF = dfF.loc[dfF.Dx != "control"]
        if caseDF.empty:
            caseDF = None
    else:
        controlDF = dfF.copy()
        caseDF = None

    # Processing

    from neuroHarmonize import harmonizationLearn, saveHarmonizationModel, harmonizationApply

    # TODO: add possible imputation and threshold for which ROIs to keep/drop

    rois2use = [c for c in roi_columns if c in controlDF.columns]
    missing_in_action = [c for c in roi_columns if c not in controlDF.columns]
    if len(missing_in_action) > 0:
        txtlog.stdO_file("\nOnly harmonizing subset of ROIs (%i) found in reference\n" % (len(roi_columns) - len(rois2use)))
        log["FA"]["missing_ROIs"] = missing_in_action

    Ncontrols = controlDF.shape[0]
    controlDFfilt = controlDF.dropna(subset=rois2use+covars)
    if controlDFfilt.shape[0] < Ncontrols:
        txtlog.stdO_file("\nDropping %i controls due to NA values\n" %(Ncontrols - controlDFfilt.shape[0]))
    
    dfRef = dfR[["Age", "Sex"] + rois2use]
    if "reference" in controlDFfilt.SITE.unique():
        raise ValueError("Cannot have a SITE named 'reference'")
    else:
        dfRef.loc[:, "SITE"] = "reference"

    dfComb = pd.concat([controlDFfilt[["subjectID"]+covars+rois2use], dfRef[covars+rois2use]], ignore_index=True)
    
    CoRarray = np.array(dfComb[rois2use])
    CoRcovars = dfComb[covars]

    run_msg_1='''
    
    ###################################
    Applying ComBat-GAM to Controls now
    ###################################

    '''

    txtlog.stdO_file(run_msg_1)
    gam_model, control_data_adj = harmonizationLearn(CoRarray, CoRcovars, smooth_terms=['Age'], ref_batch="reference")
    gam_adj_DF = pd.concat([
        dfComb[["subjectID"]].reset_index(drop=True),
        CoRcovars.reset_index(drop=True),
        pd.DataFrame(control_data_adj, columns=rois2use)],
        axis=1)
    gam_adj_DF = gam_adj_DF.loc[gam_adj_DF.SITE != "reference"]

    if caseDF is not None:
        print(caseDF.shape)
        Ncases = caseDF.shape[0]
        caseDFfilt = caseDF.dropna(subset=rois2use+covars)
        if caseDFfilt.shape[0] < Ncases:
            txtlog.stdO_file("\nDropping %i cases due to NA values\n" %(Ncases - caseDFfilt.shape[0]))
    
        # CaComb = caseDFfilt[["subjectID"]+covars+rois2use]
        CaComb = pd.concat([caseDFfilt[["subjectID"]+covars+rois2use], dfRef[covars+rois2use]], ignore_index=True)
        CaRarray = np.array(CaComb[rois2use])
        CaRcovars = CaComb[covars]
        run_msg_2='''
        
        ###################################
        Applying ComBat-GAM to Cases now
        ###################################

        '''

        txtlog.stdO_file(run_msg_2)
        case_data_adj = harmonizationApply(CaRarray, CaRcovars, gam_model)

        case_adj_DF = pd.concat([
            caseDFfilt[["subjectID"]].reset_index(drop=True),
            CaRcovars.reset_index(drop=True),
            pd.DataFrame(case_data_adj, columns=rois2use)],
            axis=1)
        case_adj_DF = case_adj_DF.loc[case_adj_DF.SITE != "reference"]
        gam_adj_DF = pd.concat([gam_adj_DF, case_adj_DF], ignore_index=True)

    dfO = dfI.set_index("subjectID")
    adj_df = gam_adj_DF.set_index("subjectID")
    dfO = dfO.loc[dfO.index.isin(adj_df.index)]
    for k, v in roi_dict.items():
        adj_df.rename(columns={v: k}, inplace=True)
        dfO[k] = adj_df[k]

    # Wrapping Up

    dfO.to_csv(os.path.join(outdir, "harmonized_FA_%s.csv" %reference), index_label="subjectID")
    saveHarmonizationModel(gam_model, os.path.join(outdir, "harmonized_FA_%s.model" %reference))
    with open(os.path.join(outdir, "harmonized_FA_%s.json" % reference), 'w') as f:
        json.dump(log, f)

    qcdir = os.path.join(outdir, "QC_images")
    if not os.path.exists(qcdir):
        os.makedirs(qcdir)
    txtlog.stdO_file("\nMaking QC Images now\n")
    efu.qc_images(rois2use, dfR, dfF, gam_adj_DF, qcdir)
    txtlog.stdO_file("\n\nSuccessful Completion!\n")
    txtlog.to_file(os.path.join(outdir, "harmonized_FA_%s.txt" %reference))

    return 

@eharmonize.command()
@click.option("--incsv", metavar="CSV",
        help="File path of the CSV file with the necessary covariates and FA measures")
@click.option("--model", metavar="MODEL", help="Model file as output by initial harmonization")
@click.option("--log", metavar="JSON", help="JSON file from initial harmonization")
@click.option("--outdir", metavar="DIR", help="Directory to write out outputs to")
@click.option("--metric", default="FA", show_default=True, metavar="METRIC", help="(optional)")
def apply_harmonization(incsv, model, outdir, log, metric):

    from neuroHarmonize import loadHarmonizationModel, harmonizationApply
    
    # Settings Check

    if not os.path.exists(incsv):
        raise OSError("--incsv input does not exist:\n%s" %incsv)

    if not os.path.exists(model):
        raise OSError("--model input does not exist:\n%s" %model)

    if not os.path.exists(log):
        raise OSError("--log input does not exist:\n%s" %log)

    with open(log, 'r') as f:
        dictI = json.load(f)
    reference = dictI["reference"]

    global txtlog
    txtlog = efu.logstr()
    log2 = txtlog.run_settings
    log2["reference"] = reference

    dfI = pd.read_csv(incsv)
    covars = ["Age", "Sex", "SITE"]
    if "SITE" not in dfI.columns:
        sites = [k for k in dictI.keys() if "SITE" in k]
        if len(sites) == 1:
            dfI.loc[:, "SITE"] = "enigma"
        else:
            raise ValueError("\nThe initial model was run with multiple sites.\nPlease include a SITE column indicating which site(s) the new data has come from.\n")

    if "Dx" in dfI.columns:
        dfI.loc[dfI.Dx.notnull(), "Dx"] = dfI.loc[dfI.Dx.notnull(), "Dx"].str.strip().str.lower()
    cohort_info = efu.input_check(dfI, covars[:-1])
    log2.update(cohort_info)

    pkg_settings, dfR = efu.load_version(reference, metric)
    roi_columns = dfR.columns[3:] 
    dfF, roi_dict, outM = efu.column_match(dfI, roi_columns, metric)
    txtlog.stdO_file(outM)
    txtlog.stdO_file("\nUsing ROIS:\n%s\n" % ", ".join(roi_dict.keys()))
    log2["FA"] = {"ROIs": list(roi_dict.keys())}

    # Subject Check

    age_excludes, outM = efu.age_check(dfF, pkg_settings)
    txtlog.stdO_file(outM)
    dfF = dfF.loc[~age_excludes]

    # Processing

    Ndata = dfF.shape[0]
    rois2use = list(roi_dict.values())
    newDFfilt = dfF.dropna(subset=rois2use+covars)
    if newDFfilt.shape[0] < Ndata:
        txtlog.stdO_file("\nDropping %i cases due to NA values\n" %(Ndata - newDFfilt.shape[0]))

    dfRef = dfR[["Age", "Sex"] + rois2use]
    if "reference" in newDFfilt.SITE.unique():
        raise ValueError("Cannot have a SITE named 'reference'")
    else:
        dfRef.loc[:, "SITE"] = "reference"

    newComb = pd.concat([newDFfilt[["subjectID"]+covars+rois2use], dfRef[covars+rois2use]], ignore_index=True)
    # newComb = newDFfilt[["subjectID"]+covars+rois2use]
    newRarray = np.array(newComb[rois2use])
    newRcovars = newComb[covars]
    run_msg='''
    
    ###################################
    Applying ComBat-GAM to all data now
    ###################################

    '''

    txtlog.stdO_file(run_msg)
    gam_model = loadHarmonizationModel(model)
    new_data_adj = harmonizationApply(newRarray, newRcovars, gam_model)

    new_adj_DF = pd.concat([
        newDFfilt[["subjectID"]].reset_index(drop=True),
        newRcovars.reset_index(drop=True),
        pd.DataFrame(new_data_adj, columns=rois2use)],
        axis=1)
    new_adj_DF = new_adj_DF.loc[new_adj_DF.SITE != "reference"]

    dfO = dfI.set_index("subjectID")
    new_adj_DF.set_index("subjectID", inplace=True)
    dfO = dfO.loc[dfO.index.isin(new_adj_DF.index)]
    for k, v in roi_dict.items():
        dfO.loc[:, k] = new_adj_DF[v]

    # Wrapping Up

    batch_msg = """
    Figuring out which batch this is
    """
    txtlog.stdO_file(batch_msg)
    batch_no = 2
    csvO = os.path.join(outdir, "batch%02d_harmonized_FA_%s.csv" %(batch_no, reference))
    while os.path.exists(csvO):
        batch_no += 1
        csvO = os.path.join(outdir, "batch%02d_harmonized_FA_%s.csv" %(batch_no, reference))
    txtlog.stdO_file("\nSettled on batch %02d\n" %batch_no)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    base_file_name = "batch%02d_harmonized_FA_%s" %(batch_no, reference)
    dfO.to_csv(os.path.join(outdir, base_file_name+".csv"), index_label="subjectID")
    with open(os.path.join(outdir, base_file_name+".json"), 'w') as f:
        json.dump(log, f)

    qcdir = os.path.join(outdir, "batch%02d_QC_images" %batch_no)
    if not os.path.exists(qcdir):
        os.makedirs(qcdir)
    txtlog.stdO_file("\nMaking QC Images now\n")
    efu.qc_images(rois2use, dfR, dfF, new_adj_DF, qcdir)
    txtlog.stdO_file("\n\nSuccessful Completion!\n")
    txtlog.to_file(os.path.join(outdir, base_file_name+".txt"))

    return

# @eharmonize.command()
# @click.option()
# def harmonize_all(indir, outdir):

if __name__ == "__main__":
    eharmonize()
