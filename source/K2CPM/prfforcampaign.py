import numpy as np

from campaigngridradec2pix import CampaignGridRaDec2Pix
from prfdata import PrfData


class PrfForCampaign(object):
    """gives PRF values for a star and set of pixels in given (sub-)campaign"""

    def __init__(self, campaign=None, grids=None, prf_data=None):
        """grids must be of CampaignGridRaDec2Pix type and prf_data must be of PrfData type"""
        if not isinstance(grids, CampaignGridRaDec2Pix):
            raise TypeError("wrong type of of grids argument ({:}), CampaignGridRaDec2Pix expected".format(type(grids)))
        if not isinstance(prf_data, PrfData):
            raise TypeError("wrong type of of prf_data argument ({:}), PrfData expected".format(type(prf_data)))
            
        if campaign is None:
            self.campaign = None
        else:
            self.campaign = int(campaign)

        self.grids = grids
        self.prf_data = prf_data

    def mean_position(self, ra, dec):
        """for given (RA,Dec) calculate position for every epoch
        and report weighted average"""
        return self.grids.mean_position(ra=ra, dec=dec)

    def apply_grids_and_prf(self, ra, dec, pixels, bjds=None, fast=True):
        """For star at given (RA,Dec) try to calculate its positions for 
        epochs (if possible) and for every epoch predict PRF value for every 
        pixel in pixels (type: list). fast=True means claculations will be 
        much faster and just slightly less accurate. 
        
        Returns:
          failed - number of epochs for which there are no grids
          mask - numpy.array of shape (len(bjds),) that gives the mask for 
            calculated epochs; you will later use bjds[mask] etc. in your code.
          out_prfs - numpy.array of shape (len(bjds), len(pixels)) that gives 
            prf values for all epochs and pixels; remember to use 
            out_prfs[mask]
        """
        if not isinstance(ra, float) or not isinstance(dec, float):
            raise TypeError('wrong types of RA,Dec in apply_grids_and_prf(): {:} and {:}; 2 floats expected'.format(type(ra), type(dec)))

        (positions_x, positions_y) = self.grids.apply_grids(ra=ra, dec=dec)

        mask = np.copy(self.grids.mask)
        out_prfs = np.zeros(shape=(len(mask), len(pixels)))
        for (i, mask_) in enumerate(mask):
            if not mask_:
                continue
            out_prfs[i] = self.prf_data.get_interpolated_prf(
                            positions_x[i], positions_y[i], pixels, fast=fast)
        
        self._positions_x = positions_x
        self._positions_y = positions_y
        
        failed = len(mask) - sum(mask)
        
        return (failed, mask, out_prfs)

    @property
    def grids_bjd_short(self):
        """return BJD array for grids (2450000 is subtracted)"""
        return self.grids.bjd - 2450000.
        
    def get_mask_and_index_for_time(self, time):
        """map the two vectors to very close values;
        returns mask to be applied to input time: time[mask]
        and index vector to be applied to apply_grids_and_prf() results: prfs[index]"""
        index = []
        mask = np.ones(len(time), dtype='bool')
        
        for (i, t) in enumerate(time):
            if t < 2450000.:
                t += 2450000.
            try:
                index.append(self.grids.index_for_bjd(t))
            except:
                mask[i] = False
                
        return (mask, index)
