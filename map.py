import numpy as np
from numpy import ma

import moby2
import moby2.scripting.products as products

from todloop import Routine

from utils import *


class PrescanTOD(Routine):
    def __init__(self, **params):
        Routine.__init__(self)
        self.params = params

    def initialize(self):
        # Expand mind of pointing model
        user_config = moby2.util.get_user_config()
        moby2.pointing.set_bulletin_A(params=user_config.get('bulletin_A_settings'))

        self.tod_names = {}   # translate from tod_files to tod names
        self.tod_info = {}    # key is tod.info.name
        self.map_lims = {}    # key is map_name

    def execute(self, store):
        tod_file = self.get_name()

        # Load tod info
        ti = moby2.scripting.get_tod_info({'filename': tod_file})

        # Start with cuts
        cuts_tod = products.get_cuts(self.params['tod_cuts'], tod=ti.basename)
        if self.params.get('position_cuts') is not None:
            cuts_pos = products.get_cuts(params['position_cuts'],
                                         tod=ti.basename)
            cuts_pos.merge_tod_cuts(cuts_tod, cut_missing=True)
        else:
            cuts_pos = cuts_tod

        # Update load_args?
        load_args = self.params.get('load_args', {})

        use_det_cuts = self.params.get('use_cuts')
        use_sample_cuts = use_det_cuts or \
                          self.params.get('use_sample_cuts')

        if use_sample_cuts:
            load_args['start'] = cuts_tod.sample_offset
            load_args['end'] = cuts_tod.sample_offset + cuts_tod.nsamps

        if use_det_cuts:
            load_args['det_uid'] = cuts_tod.get_uncut(det_uid=True)

        det_uid_args = self.params('det_uid')
        if det_uid_args is not None:
            if det_uid_args.get('filename'):
                det_uid_args['filename'] = det_uid_args['filename'].\
                                           format(**ti.get_dict())
                uid_list = moby2.detectors.DetectorList.from_dict(
                    inputm.get_filename(det_uid_args['filename']))
                load_args['det_uid'] = uid_list.det_uid
            else:
                raise

        # Load TOD pointing
        self.logger.info("Preliminary load of %s" % tod_file)
        load_args.update({'filename': tod_file,
                          'read_data': False})
        tod = moby2.scripting.get_tod(load_args)
        moby2.tod.repair_pointing(tod)

        # Trim the cuts to match the TOD.
        for cuts in [cuts_tod, cuts_pos]:
            if cuts.sample_offset > tod.info.sample_index:
                raise ValueError, "load_args.start earlier than cuts file start."
            if cuts.nsamps < tod.nsamps:
                raise ValueError, "load_args.end goes beyond end of cuts."
            # Shift.
            for c in cuts.cuts:
                c[:] -= tod.info.sample_index - cuts.sample_offset
                c.nsamps = tod.nsamps
            cuts.nsamps = tod.nsamps
            cuts.sample_offset = tod.info.sample_index
            cuts.cuts = [c.get_collapsed() for c in cuts.cuts]

        info['cuts'] = (cuts_tod, cuts_pos)

        all_det_uid = tod.info.array_data.select_outer()

        try:
            info['name'] = tod.info.name
        except:
            info['name'] = '...' + tod_file[-40:]
        self.tod_names[tod_file] = info['name']

        # Detector offsets
        info['fplane'] = products.get_focal_plane(
            self.params['pointing'], tod_info=tod.info, det_uid=all_det_uid,
            tod=tod)

        info['taus'] = products.get_time_constants(self.params['time_constants'], tod.info)

        info['cal_vect'] = products.get_calibration(
            self.params['calibration'], tod.info, det_uid=all_det_uid)

        # HWP?
        self.logger.info("Loading HWP angles...")
        hwp_params = self.params.get('hwp_angle')
        if hwp_params is not None:
            info['hwp'] = products.get_hwp_angles(hwp_params, tod=tod)
            info['hwp'] *= np.pi/180 # convert to radians.
            # Possible sign flip?
            info['hwp_signed'] = info['hwp'] * hwp_params.get('sign_convention', +1)
        else:
            self.logger.info('... no hwp angles loaded.')
            info['hwp_signed'] = None
            info['hwp'] = None

        # Get a wand for each map.
        info['wands'] = {}

        self.logger.info("Setting up projections for %s" % info['name'])
        for map_name, map_params in self.params['maps']:
            self.logger.info("... for map %s" % map_name)
            tod.fplane = info['fplane']    # Grimace

            source_name = map_params.get('source_name')
            if source_name is not None:
                self.logger.info('User specified source_name=%s' % str(source_name))
            else:
                matched_sources = moby2.ephem.get_sources_in_patch(
                    tod=tod, source_list=self.params.get(('planet_cuts', 'source_list')))
                if len(matched_sources) == 0:
                    self.logger.info('No sources found, proceeding anyway.')
                else:
                    source_name = matched_sources[0]
                    if len(matched_sources) > 1:
                        self.logger.info('Found sources %s; using %s' %
                                         (matched_sources, source_name))

            if map_params['coords'] in ['source_scan', 'source_centered']:
                ra_src, dec_src = moby2.ephem.get_source_coords(
                    source_name, tod.ctime.mean())
                wand = moby2.pointing.ArrayWand.for_tod_source_coords(
                    tod, ref_coord=(ra_src, dec_src),
                    scan_coords=map_params['coords']=='source_scan',
                    polarized=True, hwp=info.get('hwp_signed'))
                info['wands'][map_name] = wand

            elif map_params['coords'] in ['J2000', 'equatorial']:
                wand = moby2.pointing.ArrayWand.for_tod(
                    tod, hwp=info.get('hwp_signed'))
                info['wands'][map_name] = wand

            else:
                raise ValueError, "unknown map_params['coords'] = %s" % \
                    map_params['coords']

            # Also get the pointing limits for this tod, this map.
            s = tod.fplane.mask
            x0, y0 = wand.fcenter.x[0], wand.fcenter.y[0]
            x, y = tod.fplane.x[s] - x0, tod.fplane.y[s] - y0
            r_max = np.max(x**2 + y**2)**.5 * 1.05
            theta = np.arange(0., 360., 20) * np.pi/180
            cplane = moby2.pointing.FocalPlane(x=r_max*np.cos(theta)+x0,
                                               y=r_max*np.sin(theta)+y0)
            x, y = wand.get_coords(cplane)[:2]
            x, y = -x * 180./np.pi, y * 180./np.pi

            if not map_name in self.map_lims:
                self.map_lims[map_name] = [(x.min(), x.max()),
                                           (y.min(), y.max())]
            else:
                # Merge
                self.map_lims[map_name] = [ (min(a[0], _x.min()),
                                             max(a[1], _x.max()))
                                            for a, _x in zip(self.map_lims[map_name], [x,y]) ]
            del wand

        self.tod_info[info['name']] = info
        del tod
