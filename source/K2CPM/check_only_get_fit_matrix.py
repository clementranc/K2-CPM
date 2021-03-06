import sys
import numpy as np

import k2_cpm
import epic


def get_tpf_file_name(directory, epic_id, campaign):
    """get the name of the TPF file"""
    return '{:}/ktwo{:}-c{:}_lpd-targ.fits.gz'.format(directory, epic_id, campaign)

def get_epoch_mask_for_epic_id(directory, epic_id, campaign):
    """reads the file and gets epoch_mask from it"""
    file_name = get_tpf_file_name(directory, epic_id, campaign)
    #print(epic_id)
    sys.stdout.flush()
    epic.load_tpf(epic_id, campaign, directory)
    tpf = k2_cpm.Tpf(file_name)
    return tpf.epoch_mask
    
def get_predictor_epoch_mask(tpfs_epic_list_file, directory, campaign):
    """take a list of tpf files and construc predictor_epoch_mask using them"""
    predictor_epoch_mask = None
    with open(tpfs_epic_list_file) as epic_list:
        for line in epic_list.readlines():
            epoch_mask = get_epoch_mask_for_epic_id(directory, int(line), campaign)
            if predictor_epoch_mask is None:
                predictor_epoch_mask = np.ones(epoch_mask.shape, dtype=bool)
            predictor_epoch_mask &= epoch_mask
    return predictor_epoch_mask
    
def read_true_false_file(file_name):
    """reads file with values True or False"""
    parser = {'TRUE': True, 'FALSE': False}
    out = []
    with open(file_name) as in_file:
        for line in in_file.readlines():
            out.append(parser[line[:-1].upper()])
    return np.array(out)

if __name__ == '__main__':
    # no. 1 settings:    
    n_test = 1
    
    # common settings:
    l2 = 1000.
    epic_id = "200071074"
    directory = 'tpf/'
    campaign = 92 
    pixel = [883, 670]
    pre_matrix_file = "../test/output/{:}-pre_matrix.dat".format(n_test)
    tpfs_epic_list_file = "../test/output/1-epic.dat" # Content is the same for each file (1,2,3,4).
    predictor_epoch_mask_file = "1-predictor_epoch_mask.dat" # Content is the same for each file (1,2,3,4).
    fit_matrix_file_name = '{:}-fit_matrix_comp.dat'.format(n_test)
     
    pre_matrix = np.loadtxt(pre_matrix_file)

    epic.load_tpf(epic_id, campaign, directory)
    tpf_file_name = get_tpf_file_name(directory, epic_id, campaign)
    tpf = k2_cpm.Tpf(tpf_file_name)    
    
    predictor_epoch_mask = read_true_false_file(predictor_epoch_mask_file)
    
    x = pixel[0]-tpf.ref_row
    y = pixel[1]-tpf.ref_col
    index = tpf.get_index(x,y)
    
    result = k2_cpm.get_fit_matrix_ffi(
                     tpf.flux[:,index], tpf.epoch_mask, 
                     pre_matrix, 
                     predictor_epoch_mask, l2, tpf.time, 0, None)
    (target_flux, predictor_matrix, none_none, l2_vector, time, epoch_mask, data_mask) = result
    np.savetxt(fit_matrix_file_name, predictor_matrix, fmt='%.8f')
    
    