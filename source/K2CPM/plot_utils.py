import numpy as np
import matplotlib.pyplot as plt


def plot_matrix_subplots(figure, time, matrix, same_y_axis=True, data_mask=None):
    """
    Plot given 3D matrix in subpanels. Note that 3rd component of matrix.shape 
    must be the same as time.size i.e., matrix.shape[2]==time.size
    
    TO_BE_DONE: add options to specify symbols, add errorbars, add mask options
    
    example usage:
    import matplotlib.pyplot as plt
    fig = plt.gcf()
    fig.set_size_inches(50,30)
    plot_matrix(fig, time, matrix)
    plt.savefig("file_name.png")
    plt.close()
    
    to construct input data use e.g.:
    tpfdata.TpfData.directory = TPF_DIRECTORY_NAME
    tpf = tpfdata.TpfData(epic_id=SOME_epic_id, campaign=SOME_campaign)
    time = tpf.jd_short
    matrix = tpf.get_fluxes_for_square(pix_y, pix_x, half_size=3) # 3 gives 7x7 sublots
    """
    x_lim_expand = 0.06
    y_lim_expand = 0.08

    if data_mask is not None:
        y_lim = [np.nanmin(matrix[:,:,data_mask]), np.nanmax(matrix[:,:,data_mask])]
    else:    
        y_lim = [np.nanmin(matrix), np.nanmax(matrix)]
    d_y_lim = y_lim[1] - y_lim[0]  
    y_lim[0] -= d_y_lim * y_lim_expand / 2.
    y_lim[1] += d_y_lim * y_lim_expand / 2.

    (i_max, j_max, _) = matrix.shape
    panels = np.flipud(np.arange(i_max*j_max).reshape(i_max, j_max)) + 1
    if data_mask is not None:
        time = time[data_mask]
    d_time = (max(time) - min(time)) * x_lim_expand / 2.
    x_lim = [min(time)-d_time, max(time)+d_time]    

    for i in range(i_max):
        for j in range(j_max):
            ax = plt.subplot(i_max, j_max, panels[i, j])
            y_axis = matrix[i][j]
            if data_mask is not None:
                y_axis = y_axis[data_mask]
                
            ax.plot(time, y_axis, '.k')
            
            if i != 0:
                ax.get_xaxis().set_visible(False)
            if j != 0:
                ax.get_yaxis().set_visible(False)
            if same_y_axis:
                ax.set_ylim(y_lim)
            ax.set_xlim(x_lim)
    figure.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)

def construct_matrix_from_list(pixel_list, time_series_list):
    """
    take a list of pixels (pixel_list) and a list of corresponding 
    flux values (time_series_list) and make matrix that can be given to 
    plot_matrix_subplots
    """
    (x_0, y_0) = np.min(pixel_list, axis=0)
    (x_range, y_range) = np.max(pixel_list, axis=0) - np.array([x_0, y_0]) + 1

    value_length = len(time_series_list[0])
    for values in time_series_list:
        if len(values) != value_length:
            raise ValueError('construct_matrix_from_list() - all ' +
                    'time series vectors have to be of the same ' +
                    'length')

    matrix = np.zeros((x_range, y_range, value_length))

    for (i, (x, y)) in enumerate(pixel_list):
        matrix[x-x_0, y-y_0, :] = time_series_list[i]

    return matrix

