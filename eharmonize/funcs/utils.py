import pandas as pd
import pkg_resources
import numpy as np
import json
import os
import sys 
from datetime import datetime

class logstr():
    
    def __init__(self):
        now = datetime.now()
        dt_string = now.strftime("%B %d, %Y %H:%M:%S")

        self.message = "\n\nThis is a text file log for the eharmonize program.\n"
        self.run_settings = {
            "user": os.environ["USER"],
            "date": dt_string,
            "command": " ".join(sys.argv),
            }

    def stdO_file(self, message):
        self.message += message
        print(message)

    def to_file(self, outfile):
        log2write = "Runtime Information\n"
        for k, v in self.run_settings.items():
            log2write += "* %s: %s\n" %(k.title(), v)

        log2write += self.message

        with open(outfile, 'w') as f:
            f.write(log2write)

    # untested
    def abort(self, errortype, message, outfile):
      self.message += message
      self.to_file(outfile)
      raise errortype(message)

def input_check(dfI, covars):
    for covar in covars:
        if covar not in dfI.columns:
            raise KeyError("%s is a required column. Please check columns to ensure it is included/spelled correctly." %covar)

    cohort_info = {}
    sex_values = dfI.Sex.fillna("NA").value_counts()
    for k, v in sex_values.iteritems():
        if k == 0:
            cohort_info["N_Females"] = v
        elif k == 1:
            cohort_info["N_Males"] = v
        else:
            cohort_info["N_Sex_%s" % k] = v
    
    site_values = dfI.SITE.fillna("NA").value_counts()
    for k, v in site_values.iteritems():
        cohort_info["N_SITE_%s" %k] = v

    if "Dx" in dfI.columns:
        control_values = dfI.Dx.fillna("NA").value_counts()
        for k,v in control_values.iteritems():
            if k == "case":
                cohort_info["N_cases"] = v
            elif k == "control":
                cohort_info["N_controls"] = v
            else:
                cohort_info["N_Dx_%s" %k] = v
    else:
        cohort_info["N_controls"] = dfI.shape[0]

    return cohort_info

def load_version(reference, metric):
    with open(pkg_resources.resource_filename(__name__, '../data/reference_meta.json'), 'r') as f:
        pkg_meta = json.load(f)
    pkg_settings = pkg_meta[reference]

    dfR = pd.read_csv(pkg_resources.resource_filename(__name__, '../data/%s' %pkg_settings[metric]["filepath"]))

    return pkg_settings, dfR

def age_check(dfI, pkg_settings):
    outMessage = "\nChecking age range of input data\n"
    age_min = pkg_settings["FA"]["age_min"]
    too_young = dfI.Age < age_min
    if too_young.sum() > 0:
        outMessage += "\nFound %i subjects younger than %i years old. Excluding them.\n" %(too_young.sum(), age_min)

    age_max = pkg_settings["FA"]["age_max"]
    too_old = dfI.Age > age_max
    if too_old.sum() > 0:
        outMessage += "\nFound %i subjects older than %i years old. Excluding them.\n" %(too_old.sum(), age_max)

    age_excludes = (too_young | too_old)
    return age_excludes, outMessage

def column_match(dfI, roi_columns, metric):
    rois = roi_columns.str.split("_%s" %metric).str[0]
    roi_dict = {}
    tapetum = False
    outMessage = "\nChecking ROI names of input data\n"

    for c in dfI.columns:
        if c.startswith("TAP"):
            tapetum = True
        if c in ["AverageMD", "AverageAD", "AverageRD"]:
            new_name = c.replace("Average", "AverageFA_")
        else:
            new_name = c.replace("-", ".").replace("/", ".")
        if new_name in roi_columns:
            roi_dict[c] = new_name
        elif new_name in rois:
            roi_dict[c] = new_name + "_%s" % metric

    if tapetum:
        # TODO: actually switch names away from ENIGMA to JHU
        outMessage += "\nDetected Tapetum columns, which are not included in the ENIGMA-DTI version\n"
        unc_cols = dfI.columns[dfI.columns.str.startswith("UNC")]
        unc_replace = unc_cols.str.replace("UNC", "IFO")
        if len(unc_cols) > 0:
            outMessage += "Renaming UNC columns with IFO to match ENIGMA reference\n"
            for c, r in zip(unc_cols, unc_replace):
                new_name = r.replace("-", ".")
                if new_name in roi_columns:
                    roi_dict[c] = new_name
                elif new_name in rois:
                    roi_dict[c] = new_name + "_%s" % metric
        tap_cols = dfI.columns[dfI.columns.str.startswith("TAP")]
        tap_replace = tap_cols.str.replace("TAP", "UNC")
        if len(tap_cols) > 0:
            outMessage += "Renaming TAP columns with UNC to match ENIGMA reference\n"
            for c, r in zip(tap_cols, tap_replace):
                new_name = r.replace("-", ".")
                if new_name in roi_columns:
                    roi_dict[c] = new_name
                elif new_name in rois:
                    roi_dict[c] = new_name + "_%s" % metric

    dfi = dfI.rename(columns=roi_dict)

    return dfi, roi_dict, outMessage

def qc_images(rois2use, dfR, dfI, dfO, qcdir):
    import matplotlib.pyplot as plt
    import matplotlib
    import seaborn as sns

    matplotlib.use('Agg')
    font = {'size'   : 16}
    #       'family' : 'normal',
    #       'weight' : 'bold',

    matplotlib.rc('font', **font)

    age_min = dfI.Age.min()
    age_max = dfI.Age.max()
    dfR_filt = dfR.loc[(dfR.Age < age_max) & (dfR.Age > age_min)]
    dfI.loc[:, "Data"] = "Raw"
    dfO.loc[:, "Data"] = "Harmonized"
    plot_data = pd.concat([dfI, dfO], ignore_index=True)
    plot_data.dropna(subset=["Age"], inplace=True)
    plot_data.loc[:, "Age_int"] = plot_data["Age"].astype(int)

    # TODO: add case/control distinction
    # TODO: add ENIGMA LUT so full ROI names can be spelled out
    # lutfile = pkg_resources.resource_filename(__name__, '../data/ENIGMA_LUT.txt')
    # lut = pd.read_csv(lutfile, delimiter="\t")

    for roi in rois2use:
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # TODO: switch to regplot with lowess=True
        # and add scatterplot
        sns.lineplot(data=dfR, x="Age", y=roi,
                alpha = 0.5, ax=ax, color="black")
        sns.lineplot(data=plot_data, x="Age_int", y = roi,
                hue="SITE", style="Data", ax=ax)
        # sns.scatterplot(plot_data, x="Age", y=roi, 
        #         hue="SITE",style="Data", alpha=0.5, ax=ax)
        # fig, ax = plt.subplots(1, 2 figsize=(10, 5))    
        # ax[0].plot(dfR_filt.Age, dfR_filt[roi],
        #         'k.', alpha = 0.5)
        # ax[0].set_title="ROI vs Age (years)"
        # ax[1].hist(dfR_filt[roi], bins=50, alpha=0.5)
        # sns.histplot(plot_data, x=roi, bins=10, hue="SITE",
        #         ax=ax[1])
        ax.set_title("ROI vs Age (years)")
        ax.set_xlabel("Age (years)")
        ax.legend(bbox_to_anchor=(1.1, 1.05))
        plt.tight_layout()
        plt.savefig(os.path.join(qcdir, "%s.png" %roi))
    return
