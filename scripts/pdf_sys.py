#!/usr/bin/env python
import collections
import subprocess
import numpy as np

import h5py
import lhapdf

from lhapdf_lookup import pset_name

variation_keyname = 'variations'
binid_keyname = 'binids'
h1d_keyname = 'Histo1D'
h2d_keyname = 'Histo2D'
p1d_keyname = 'Profile1D'
s1d_keyname = 'Scatter1D'
s2d_keyname = 'Scatter2D'
cnt_keyname = 'Counter'
all_data_keynames = [h1d_keyname, h2d_keyname, p1d_keyname, s1d_keyname, s2d_keyname]
all_hist_keynames = [h1d_keyname, h2d_keyname, p1d_keyname]
sumw_keyname = 'sumw'
yval_keyname = 'yval'
invalid_number = -999999

def install_pdfs(pdf_list):
    # installed_pdfs = subprocess.check_output(['lhapdf', 'list', '--installed']).decode('utf-8')
    installed_pdfs = ""
    # print("Installed PDFs", installed_pdfs)
    installed_pdfs = installed_pdfs.split('\n')
    for pdf in pdf_list:
        if pdf not in installed_pdfs:
            print("Installing:", pdf)
            subprocess.call(['lhapdf', 'install', pdf])
        else:
            # print(pdf, "already installed")
            pass


def get_pdf_id(variation):
    return int(variation.split('_')[-1][3:])


def get_obs_values(input_h5):
    n_variations = input_h5[variation_keyname].shape[0]
    n_bins = input_h5[binid_keyname].shape[0]
    obs_values = np.full((n_bins, n_variations), fill_value=invalid_number, dtype=np.float32)
    sumw = input_h5[sumw_keyname]
    yval = input_h5[yval_keyname]

    iline = 0
    for hist_keyname in all_hist_keynames:
        for hist in input_h5[hist_keyname]:
            bin_idx = input_h5[hist_keyname+'/'+hist][:]
            obs_values[bin_idx] = sumw[bin_idx]
            if iline == 0:
                print(hist)
                print(bin_idx, len(bin_idx))
                print(sumw[bin_idx])
                print(sumw[bin_idx].shape)
            iline += 1
    
    for s1 in input_h5[s1d_keyname]:
        bin_idx = input_h5[s1d_keyname+'/'+s1][:]
        obs_values[bin_idx, :] = yval[bin_idx, :]

    print(obs_values[1])
    return obs_values


def get_sys(input_h5_name):
    input_h5 = h5py.File(input_h5_name, 'r')

    print(list(input_h5.keys()))
    print(list(input_h5['Histo1D'].keys()))
    # print(input_h5[h1d_keyname+'/|MC_GENERIC|Pt'][:])
    # for key in input_h5.keys():
    #     print(key)

    observables = []
    for key_name in all_data_keynames:
        values = input_h5[key_name]
        for value in values:
            observables.append(value)
    print("total {} observables".format(len(observables)))

    variations = input_h5[variation_keyname]
    n_variations = len(variations)
    print("total {} variations".format(n_variations))

    binids = input_h5[binid_keyname]
    print(binids.shape)
    print(binids[0])
    n_bins = len(binids)
    print("total {} bins".format(n_bins))


    nominal_pdf_name = "none"
    nominal_pdf_id = -1
    n_exp_internal_pdf_vars = -1
    pdfs_variations =  collections.defaultdict(lambda: 0)
    pdfs_variation_idx =  collections.defaultdict(list)
    external_pdfs = []
    for idx, variation in enumerate(variations):
        if "MUR1_MUF0.5" in variation:
            nominal_pdf_id = get_pdf_id(variation)
            nominal_pdf_name = pset_name(nominal_pdf_id)
            if not nominal_pdf_name:
                raise ValueError("Did you find nominal PDF:", nominal_pdf_id)
        if "MUR1_MUF1" in variation:
            pdf_id = get_pdf_id(variation)
            pdf_name = pset_name(pdf_id)
            pdfs_variations[pdf_name] += 1
            pdfs_variation_idx[pdf_name].append(idx)

    print("nominal PDF:", nominal_pdf_name)
    n_external_pdfs = 0
    for key,value in pdfs_variations.items():
        print("{}: {} variations".format(key, value))
        print("\t", pdfs_variation_idx[key])
        if key != nominal_pdf_name:
            n_external_pdfs += 1
    print("{} external PDFs".format(n_external_pdfs))

    install_pdfs(list(pdfs_variations.keys()))
    pset = lhapdf.getPDFSet(nominal_pdf_name)
    pdfs_ = pset.mkPDFs()
    n_exp_internal_pdf_vars = len(pdfs_)

    # evaluate the PDF sys
    obs_values = get_obs_values(input_h5)
    pdf_sys = np.full((n_bins, 3), fill_value=invalid_number, dtype=np.float32)
    nominal_pdf = lhapdf.getPDFSet(nominal_pdf_name)
    cl = 68
    internal_pdf_sys = np.apply_along_axis(
        lambda x: nominal_pdf.uncertainty(x[pdfs_variation_idx[nominal_pdf_name]], cl).errsymm,
        1, obs_values)
    print(internal_pdf_sys[1])
    print(internal_pdf_sys.shape)
    # external PDF sets
    external_pdf_values = np.zeros((n_bins, n_external_pdfs+1), dtype=np.float32)
    inorm = -1
    idx = 0
    for key,value in pdfs_variations.items():
        print(key, value)
        if key == nominal_pdf_name:
            inorm = idx
        if value > 1:
            this_pdf = lhapdf.getPDFSet(key)
            external_pdf_values[:, idx] = np.apply_along_axis(
                lambda x: this_pdf.uncertainty(x[pdfs_variation_idx[key]], cl).central,
                1, obs_values)
        else:
            external_pdf_values[:, idx] = obs_values[:, pdfs_variation_idx[key][0]]
        idx += 1

    external_diff = external_pdf_values - external_pdf_values[:, inorm].reshape(-1, 1)
    print(np.any(np.nonzero(external_diff[:, inorm])))
    external_pdf_sys = np.amax(external_diff, axis=1)
    pdf_sys[:, 0] = internal_pdf_sys
    pdf_sys[:, 1] = external_pdf_sys
    pdf_sys[:, 2] = np.amax(pdf_sys[:, 0:2], axis=1)

    out_file = h5py.File("Rivet_pdfsys.h5", 'w')
    out_file.create_dataset("pdfsys", data=pdf_sys, dtype=np.float32)
    out_file.close()

    input_h5.close()

if __name__ == "__main__":
    get_sys('Rivet.h5')