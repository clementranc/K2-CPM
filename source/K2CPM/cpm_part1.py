from __future__ import print_function

import sys
import numpy as np
from sklearn.decomposition import PCA
import argparse

from K2CPM import matrix_xy, multipletpf, tpfdata
from K2CPM import wcsfromtpf


def run_cpm_part1(channel, campaign, num_predictor, num_pca, dis, excl, 
                    flux_lim, input_dir, pixel_list=None, train_lim=None, 
                    output_file=None, output_file_mask=None,  
                    return_predictor_epoch_masks=False, campaign_2=None):
    """
    CPM part 1 - prepare predictor matrix

    campaign_2 - if specified, then additional predictor matrix is returned 
    that uses exactly the same pixels. 
    """
# REMOVED: l2, output_dir, target_epic_num
# ADDED: output_file, output_file_mask, return_predictor_epoch_masks, channel, campaign_2
#def run_cpm_part1(target_epic_num, campaign, num_predictor, l2, num_pca, dis, excl, flux_lim, input_dir, output_dir, pixel_list=None, train_lim=None):

    #print(channel, campaign, num_predictor, num_pca, dis, excl, 
    #                flux_lim, input_dir, pixel_list, train_lim, 
    #                output_file, output_file_mask,  
    #                return_predictor_epoch_masks)

    if channel > 200000000:
        raise ValueError("Oppps... it seems you're trying to use previous " + 
                "version of run_cpm_part1(). Currently, the first parameter " +
                "is channel (i.e. 30, 31, 32 etc.), not target_EPIC_id as it " +
                "was before")
    if pixel_list is None:
        raise ValueError("Empty pixel_list in run_cpm_part1()")
    # CHANGE the part below (i.e. implement what is needed)
    if pixel_list.shape[0] != 1 and (output_file is not None or output_file_mask is not None):
        raise ValueError('\n\nCurrently we can deal with only a single pixel at a time if the output file is specified')

    if train_lim is not None:
        raise ValueError("run_cpm_part1() parameter train_lim has to be None."
            + " We keep it because of older code")
    tpfdata.TpfData.directory = input_dir

    flux_lim_step_down = 0.1
    flux_lim_step_up = 0.1
    min_flux_lim = 0.1
    n_use = 15 # THIS HAS TO BE CHANGED. !!!

    wcs = wcsfromtpf.WcsFromTpf(channel, campaign)
    m_tpfs = multipletpf.MultipleTpf(campaign=campaign)
    if campaign_2 is not None:
        m_tpfs_2 = multipletpf.MultipleTpf(campaign=campaign_2)
    else:
        m_tpfs_2 = None

    out_predictor_matrixes = []
    out_predictor_matrixes_2 = []
    out_predictor_epoch_masks = [] 
    out_predictor_epoch_masks_2 = []

    for pixel in pixel_list:
        print(pixel[0], pixel[1])
        epic_id = wcs.epic_for_pixel(row=pixel[0], column=pixel[1])
        tpf_data = m_tpfs.tpf_for_epic_id(epic_id)
        if not tpf_data.check_pixel_covered(row=pixel[0], column=pixel[1]):
            msg = 'something went wrong in CPM part 1: {:} {:} {:}'
            raise ValueError(msg.format(pixel[0], pixel[1], epic_id))
        
        (epics_to_use_all, _, _) = wcs.get_epics_around_pixel(row=pixel[0], column=pixel[1])
        epics_to_use = epics_to_use_all[:n_use] # THIS HAS TO BE CORRECTED
        m_tpfs.add_tpf_data_from_epic_list(epics_to_use, campaign)
        if campaign_2 is not None:
            m_tpfs_2.add_tpf_data_from_epic_list(epics_to_use, campaign_2)
        
        result_get_predictor = tpf_data.get_predictor_matrix(pixel[0], pixel[1], 
                                    num_predictor, dis=dis, excl=excl, 
                                    flux_lim=flux_lim, multiple_tpfs=m_tpfs, 
                                    tpfs_epics=epics_to_use, 
                                    multiple_tpfs_2=m_tpfs_2)
        if campaign_2 is not None:
            (predictor_matrix, predictor_matrix_2) = result_get_predictor
        else:
            predictor_matrix = result_get_predictor
        
        while predictor_matrix.shape[1] < num_predictor:
            n_use += 3
            epics_to_use = epics_to_use_all[:n_use]
            m_tpfs.add_tpf_data_from_epic_list(epics_to_use[-3:], campaign)
            if campaign_2 is not None:
                m_tpfs_2.add_tpf_data_from_epic_list(epics_to_use[-3:], campaign_2)
# Re-code stuff below ???
                    #low_lim = flux_lim[0] - flux_lim_step_down
                    #up_lim = flux_lim[1] + flux_lim_step_up
                    #while predictor_matrix.shape[1] < num_predictor:
                    #    old_num = predictor_matrix.shape[1]
                    #    predictor_matrix = tpf_data.get_predictor_matrix(pixel[0], pixel[1], num_predictor, dis=dis, excl=excl,
                    #                                    flux_lim=(low_lim,up_lim), 
                    #                                    multiple_tpfs=m_tpfs, tpfs_epics=epics_to_use)
                    #    low_lim = np.max(low_lim-flux_lim_step_down, min_flux_lim)
                    #    up_lim = up_lim + flux_lim_step_up
                    #    difference = predictor_matrix.shape[1] - old_num
                    #    if difference == 0:
                    #        print('no more pixel at all')
                    #        break
                    #break
            result_get_predictor = tpf_data.get_predictor_matrix(pixel[0], pixel[1], 
                                    num_predictor, dis=dis, excl=excl,
                                    flux_lim=(low_lim, up_lim), 
                                    multiple_tpfs=m_tpfs, 
                                    tpfs_epics=epics_to_use, 
                                    multiple_tpfs_2=m_tpfs_2)
            if campaign_2 is not None:
                (predictor_matrix, predictor_matrix_2) = result_get_predictor
            else:
                predictor_matrix = result_get_predictor

        if num_pca>0:
            if campaign_2 is not None:
                raise NotImplementedError("This feature has not been implemented yet")
                # I'm not sure how exactly combine 2 predictor matrixes and PCA
            pca = PCA(n_components=num_pca, svd_solver='full')
            pca.fit(predictor_matrix)
            predictor_matrix = pca.transform(predictor_matrix)
                
        out_predictor_matrixes.append(predictor_matrix)
        out_predictor_epoch_masks.append(m_tpfs.predictor_epoch_mask)
        if campaign_2 is not None:
            out_predictor_matrixes_2.append(predictor_matrix_2)
            out_predictor_epoch_masks_2.append(m_tpfs_2.predictor_epoch_mask)

        if output_file is not None: 
            if campaign_2 is not None:
                raise NotImplementedError("This feature has not been implemented yet")
            matrix_xy.save_matrix_xy(predictor_matrix, output_file)
        if output_file_mask is not None:
            if campaign_2 is not None:
                raise NotImplementedError("This feature has not been implemented yet")
            np.savetxt(output_file_mask, m_tpfs.predictor_epoch_mask, fmt='%r')

    if return_predictor_epoch_masks:
        if campaign_2 is not None:
            return (out_predictor_matrixes, out_predictor_matrixes_2, out_predictor_epoch_masks, out_predictor_epoch_masks_2)
        else:
            return (out_predictor_matrixes, out_predictor_epoch_masks)
    else:
        if campaign_2 is not None:
            return (out_predictor_matrixes, out_predictor_matrixes_2)
        else:
            return out_predictor_matrixes


def main():
    parser = argparse.ArgumentParser(description='k2 CPM')
    parser.add_argument('epic', nargs=1, type=int, help="epic number")
    parser.add_argument('campaign', nargs=1, type=int, help="campaign number, 91 for phase a, 92 for phase b")
    parser.add_argument('n_predictor', nargs=1, type=int, help="number of predictor pixels")
    parser.add_argument('n_pca', nargs=1, type=int, help="number of the PCA components to use, if 0, no PCA")
    parser.add_argument('distance', nargs=1, type=int, help="distance between target pixel and predictor pixels")
    parser.add_argument('exclusion', nargs=1, type=int, help="how many rows and columns that are excluded around the target pixel")
    parser.add_argument('input_dir', nargs=1, help="directory to store the output file")
    parser.add_argument('output_file', nargs=1, help="output predicotr matrix file")
    parser.add_argument('-p', '--pixel', metavar='pixel_list', help="path to the pixel list file that specify list of pixels to be modelled." 
                                                                    "If not provided, the whole target pixel file will be modelled")
    parser.add_argument('-t', '--train', nargs=2, metavar='train_lim', help="lower and upper limit defining the training data set")

    args = parser.parse_args()
    print("epic number: {0}".format(args.epic[0]))
    print("campaign: {0}".format(args.campaign[0]))
    print("number of predictors: {0}".format(args.n_predictor[0]))
    print("number of PCA components: {0}".format(args.n_pca[0]))
    print("distance: {0}".format(args.distance[0]))
    print("exclusion: {0}".format(args.exclusion[0]))
    print("directory of TPFs: {0}".format(args.input_dir[0]))
    print("output predictor_matrix file: {0}".format(args.output_file[0]))
    # Variables flux_lim, pixel_list, amd train_lim and used later but not printed here.

    if args.pixel is not None:
        pixel_list = np.loadtxt(args.pixel, dtype=int, ndmin=2)
        print("pixel list: {0}".format(args.pixel))
    else:
        pixel_list = None
        print("full image")
    flux_lim = (0.2, 1.5)

    if args.train is not None:
        train_lim = (float(args.train[0]), float(args.train[1]))
        print("train limit: {0}".format(train_lim)) 
    else:
        train_lim = None
        print("all data used")
        
    run_cpm_part1(args.epic[0], args.campaign[0], args.n_predictor[0], 
                    args.n_pca[0], args.distance[0], args.exclusion[0], 
                    flux_lim, args.input_dir[0], pixel_list, train_lim,
                    output_file=args.output_file[0]) 

if __name__ == '__main__':
    main()
