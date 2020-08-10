#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019 B. Biffi
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3.


import ConfigParser
import argparse
import glob
import os
import subprocess
import sys

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--img', required=True, help='Input image')
parser.add_argument('-o', '--out', required=True, help='Result segmentation')
parser.add_argument('-c', '--config', required=True, help='Config file')
parser.add_argument('--gcc', help='Use GCC ranking',
                    action='store_true')
parser.add_argument('--res', help='Keep intermediate results',
                    action='store_true')
parser.add_argument('--iso', help='Isotropic image resampling',
                    action='store_true')
args = parser.parse_args()

if args.iso:
    query_unres = args.img
    query = query_unres.replace('.nii.gz', '_resampled.nii.gz')
else:
    query = args.img

segmentation = args.out

# Read configuration file and setup search queries
config = ConfigParser.ConfigParser()
config.read(args.config)
atlasdir = config.get('atlas', 'dir')
atlaspattern = atlasdir + '/' + config.get('atlas', 'pattern')
labelpattern = atlasdir + '/' + config.get('atlas', 'labels')
maskpattern = atlasdir + '/' + config.get('atlas', 'mask')

# Basic configuration
out_folder = os.path.dirname(segmentation)  # output dir
base_name = os.path.basename(query)  # base filename
base_name = base_name.replace('.gz', '')
base_name = base_name.replace('.nii', '')
base_name = base_name.replace('.hdr', '')
base_name = base_name.replace('.img', '')

# Output tmp files
query_mask_in_final = out_folder + '/fusionedmaskinput_final_query_' + base_name + '.nii.gz'
query_mask_out      = out_folder + '/querymask_' + base_name + '.nii.gz'
query_mask_out_two  = out_folder + '/querymask_two_' + base_name + '.nii.gz'
query_t1_in_final   = out_folder + '/fusionedT1input_final_query_' + base_name + '.nii.gz'

# Query atlas list
atlas_list = glob.glob(atlaspattern)
label_list = glob.glob(labelpattern)
mask_list = glob.glob(maskpattern)
atlas_list.sort()
#print(atlas_list)
label_list.sort()
#print(label_list)
mask_list.sort()
#print(mask_list)

# Resample original image to have isotropic voxel size
if args.iso:
    if os.path.exists(query):
        print 'Image already resampled, skipping first task'
    else:
        cmd_isotropic = 'c3d ' + query_unres + ' -resample-mm 0.6x0.6x0.6mm -o ' + query
        p = subprocess.Popen(cmd_isotropic, shell=True)
        p.communicate()

# Rigid registration for rough mask
if os.path.exists(query_mask_out):
    print('Found mask skipping first task')
else:
    counter = 0
    all_masks = list()
    for an_atlas, a_mask in zip(atlas_list, mask_list):
        affine_matrix = out_folder + '/affine_ref_' + str(counter) + '_flt_' + base_name + '.txt'
        inv_affine = out_folder + '/invaffine_ref_' + str(counter) + '_flt_' + base_name + '.txt'
        resampled_mask = out_folder + '/resampled_query_' + base_name + '_to_mask_' + str(counter) + '_ini.nii.gz'
        all_masks.append(resampled_mask)

        if not os.path.exists(inv_affine):
            cmd_aladin = 'reg_aladin -ref ' + an_atlas + ' -flo ' + query + ' -rmask ' + a_mask + ' -aff ' + \
                         affine_matrix + ' -rigOnly -noSym'
            cmd_transform = 'reg_transform -ref ' + an_atlas + ' -invAffine ' + affine_matrix + ' ' + inv_affine
            print(cmd_aladin)

            p = subprocess.Popen(cmd_aladin, shell=True)
            p.communicate()

            p = subprocess.Popen(cmd_transform, shell=True)
            p.communicate()

        if not os.path.exists(resampled_mask):
            cmd_resample = ('reg_resample -ref ' + query + ' -flo ' + a_mask +
                            ' -trans ' + inv_affine + ' -res ' + resampled_mask + ' -inter 0')

            p = subprocess.Popen(cmd_resample, shell=True)
            p.communicate()

        counter += 1

    # Merge masks
    cmd_merge = ('seg_maths ' + all_masks.pop(0) + ' -merge ' + str(len(all_masks)) + ' 4 ' + ' '.join(
        ['%s' % m for m in all_masks]) + ' -tmean -thr 0.3 -bin -dil 9 ' + query_mask_out)

    p = subprocess.Popen(cmd_merge, shell=True)
    p.communicate()

    if os.path.exists(query_mask_out):
        if not args.res:
            tmp_list = glob.glob(out_folder + '/resampled_*')
            for a_tmp in tmp_list:
                os.remove(a_tmp)
    else:
        print('Something went wrong creating the first mask. Program will exit. Please check')
        sys.exit('Mask was not generated')


# Affine registration
if os.path.exists(query_mask_out_two):
    print('Affine mask found. Do we really need to do this? Skipping this bit')
else:
    counter = 0
    all_masks = list()
    for an_atlas, a_mask in zip(atlas_list, mask_list):
        inv_affine = (out_folder + '/invaffine_ref_' + str(counter) +
                      '_flt_' + base_name + '_b.txt')
        resampled_mask = (out_folder + '/resampled_mask_' + str(counter) +
                          '_mid_to_query_' + base_name + '.nii.gz')
        all_masks.append(resampled_mask)

        if not os.path.exists(inv_affine):
            cmd_aladin = ('reg_aladin -flo ' + an_atlas + ' -fmask ' + a_mask +
                          ' -ref ' + query + ' -rmask ' + query_mask_out +
                          ' -aff ' + inv_affine + ' -maxit 8')
            p = subprocess.Popen(cmd_aladin, shell=True)
            p.communicate()

        if not os.path.exists(resampled_mask):
            cmd_resample = ('reg_resample -ref ' + query + ' -flo ' + a_mask +
                            ' -trans ' + inv_affine + ' -res ' + resampled_mask +
                            ' -inter 0')
            p = subprocess.Popen(cmd_resample, shell=True)
            p.communicate()

        counter += 1

    cmd_merge = ('seg_maths ' + all_masks.pop(0) + ' -merge ' + str(len(all_masks)) +
                 ' 4 ' + ' '.join(['%s' % m for m in all_masks]) +
                 ' -tmean -thr 0.7 -bin -dil 8 ' + query_mask_out_two)

    p = subprocess.Popen(cmd_merge, shell=True)
    p.communicate()

    if os.path.exists(query_mask_out_two):
        if not args.res:
            tmp_list = glob.glob(out_folder + '/resampled_*')
            for a_tmp in tmp_list:
                os.remove(a_tmp)
    else:
        print('Something went wrong creating the second mask. \
                Program will exit. Please check')
        sys.exit('Mask was not generated')

# Non rigid registration
if os.path.exists(query_mask_in_final) and os.path.exists(query_t1_in_final):
    print('This is weird. Both 4D images exist. Performing only the lab fusion')
else:
    counter = 0
    all_final_masks = list()
    all_final_t1 = list()
    for an_atlas, a_label in zip(atlas_list, label_list):
        control_points = (out_folder + '/cpp_query_' + base_name + '_to_' +
                          str(counter) + '.nii.gz')
        inv_affine = (out_folder + '/invaffine_ref_' + str(counter) + '_flt_' +
                      base_name + '_b.txt')
        final_resampled = (out_folder + '/resampled_mask_' + str(counter) +
                           '_final_to_query_' + base_name + '.nii.gz')
        final_resampled_t1 = (out_folder + '/resampled_T1_' + str(counter) +
                              '_final_to_query_' + base_name + '.nii.gz')
        all_final_masks.append(final_resampled)
        all_final_t1.append(final_resampled_t1)

        if not os.path.exists(control_points):
            cmd_f3d = ('reg_f3d -ref ' + query + ' -rmask ' + query_mask_out_two +
                       ' -flo ' + an_atlas + ' -cpp ' + control_points + ' -aff ' +
                       inv_affine)
            p = subprocess.Popen(cmd_f3d, shell=True)
            p.communicate()

        if not os.path.exists(final_resampled_t1):
            cmd_resample = ('reg_resample -ref ' + query + ' -flo ' + an_atlas +
                            ' -trans ' + control_points + ' -res ' +
                            final_resampled_t1)
            p = subprocess.Popen(cmd_resample, shell=True)
            p.communicate()

        if not os.path.exists(final_resampled):
            cmd_resample = ('reg_resample -ref ' + query + ' -flo ' + a_label +
                            ' -trans ' + control_points + ' -res ' + final_resampled +
                            ' -inter 0')
            p = subprocess.Popen(cmd_resample, shell=True)
            p.communicate()

        counter += 1

    if args.gcc:
        num_masks = len(atlas_list) / 2
        head_mask = all_final_masks.pop(0)
        tmp_mask = out_folder + '/tempmask_' + base_name + '.nii.gz'

        cmd_merge = ('seg_maths ' + head_mask + ' -merge ' + str(len(all_final_masks)) +
                     ' 4 ' + ' '.join(['%s' % m for m in all_final_masks]) +
                     ' -bin -tmean -thr 0.2 -bin -dil 5 ' + tmp_mask)
        cmd_calc = ('seg_CalcTopNCC -target ' + query + ' -templates ' + str(len(atlas_list)) +
                    ' ' + ' '.join(['%s' % m for m in all_final_t1]) + ' -n ' +
                    str(num_masks) + ' -mask ' + tmp_mask)

        p = subprocess.Popen(cmd_merge, shell=True)
        p.communicate()

        all_final_masks.insert(0, head_mask)

        top_txt = subprocess.check_output(cmd_calc, shell=True)
        top_imgs = top_txt.split()
        top_masks = list()
        for an_img in top_imgs:
            an_index = all_final_t1.index(an_img)
            top_masks.append(all_final_masks[an_index])

        top_masks.sort()
        top_imgs.sort()

        cmd_merge = ('seg_maths ' + top_masks.pop(0) + ' -merge ' + str(len(top_masks)) +
                     ' 4 ' + ' '.join(['%s' % m for m in top_masks]) + ' ' +
                     query_mask_in_final)
        cmd_merge_t1 = ('seg_maths ' + top_imgs.pop(0) + ' -merge ' + str(len(top_imgs)) +
                        ' 4 ' + ' '.join(['%s' % m for m in top_imgs]) + ' ' +
                        query_t1_in_final)

        if not args.res:
            os.remove(tmp_mask)
    else:
        cmd_merge = ('seg_maths ' + all_final_masks.pop(0) + ' -merge ' + str(len(all_final_masks)) +
                     ' 4 ' + ' '.join(['%s' % m for m in all_final_masks]) + ' ' +
                     query_mask_in_final)
        cmd_merge_t1 = ('seg_maths ' + all_final_t1.pop(0) + ' -merge ' + str(len(all_final_t1)) +
                        ' 4 ' + ' '.join(['%s' % m for m in all_final_t1]) + ' ' +
                        query_t1_in_final)

    p = subprocess.Popen(cmd_merge, shell=True)
    p.communicate()

    p = subprocess.Popen(cmd_merge_t1, shell=True)
    p.communicate()

    if os.path.exists(query_mask_in_final) and os.path.exists(query_t1_in_final):
        if not args.res:
            tmp_list = glob.glob(out_folder + '/resampled_*')
            tmp_list.extend(glob.glob(out_folder + '/cpp*'))
            for a_tmp in tmp_list:
                os.remove(a_tmp)
    else:
        print('Something went wrong creating creating 4D images. \
                Program will exit. Please check')
        sys.exit('4D images not generated')

num_templates = len(atlas_list) / 3

cmd_fusion = ('seg_LabFusion -in ' + query_mask_in_final + ' -STEPS 6 ' + str(num_templates) + ' ' +
              query + ' ' + query_t1_in_final + ' -prop_update -out ' + segmentation +
              ' -v 1 -MRF_beta 0.8')
p = subprocess.Popen(cmd_fusion, shell=True)
p.communicate()

if not args.res and os.path.exists(segmentation):
    tmp_list = glob.glob(out_folder + '/*.txt')
    for a_tmp in tmp_list:
        os.remove(a_tmp)
    os.remove(query_t1_in_final)
    os.remove(query_mask_in_final)
    os.remove(query_mask_out_two)
    os.remove(query_mask_out)

# Keep only the largest connected component
cmd_lconcomp = 'seg_maths ' + segmentation + ' -lconcomp ' + segmentation
p = subprocess.Popen(cmd_lconcomp, shell=True)
p.communicate()
