from todloop import Routine
import cPickle, h5py, os
import numpy as np
import copy
from moby2.scripting import  pathologies_tools


class PathologyReport(Routine):
    def __init__(self, **params):
        """This routine aims to generate a Pathology Report object
        as is done for the moby2 script"""
        Routine.__init__(self)
        self.cutparam = params.get("cutparam", None)
        self.report = None

    def initialize(self):
        self.report = pathologies_tools.reportPathologies(self.cutparam)
        # update the depot file name with rank if needed
        # if mpi, label depot file with the rank to avoid racing condition
        if self.get_comm():
            self.report.depot_file += '.%d' % self.get_rank()

    def execute(self, store):
        # get obs name
        obs = self.get_name()
        # append results
        self.report.appendResult(obs)


class Summarize(Routine):
    def __init__(self, **params):
        """This routine aims to summarize the parameters derived
        from all the other routines"""
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)

    def execute(self, store):
        # initialize an empty dictionary to store results
        results = {}

        # retrieve all the calculated results
        for key in self.inputs.get('features'):
            results.update(store.get(key))

        self.logger.info("Successfully processed: %s" % results.keys())
        store.set(self.outputs.get('report'), results)


class PrepareDataLabel(Routine):
    def __init__(self, **params):
        """Prepare an HDF5 data set that contains all relevant metedata
        and detector timeseries as a basis for various ML studies
        """
        Routine.__init__(self)
        self.inputs = params.get("inputs", None)
        self._pickle_file = params.get('pickle_file', None)
        self._output_file = params.get('output_file', None)
        self._group_name = params.get('group', None)
        self._downsample = params.get('downsample', 1)
        self._remove_mean = params.get('remove_mean', False)

    def initialize(self):
        # load pickle file
        self.logger.info("Loading %s..." % self._pickle_file)
        with open(self._pickle_file, "r") as f:
            self._pickle_data = cPickle.load(f)

        # create output h5 file if it doesn't exist
        if os.path.isfile(self._output_file):
            # file exists
            self.logger.info("Found %s, updating instead" % self._output_file)
            self._hf = h5py.File(self._output_file, 'a')
        else:
            # file doesn't exist
            self.logger.info("Creating %s..." % self._output_file)
            self._hf = h5py.File(self._output_file, 'w')

        try:
            self.logger.info("Creating group %s..." % self._group_name)
            self._group = self._hf.create_group(self._group_name)
        except ValueError:
            self.logger.info("%s exists, update instead" % self._group_name)
            self._group = self._hf[self._group_name]

    def execute(self, store):
        # retrieve tod
        tod = store.get(self.inputs.get('tod'))

        # retrieve the calculated statustics
        report = store.get(self.inputs.get('report'))
        keys = report.keys()

        # remove mean if needed
        if self._remove_mean:
            report = copy.deepcopy(report)
            self.logger.info("Remove means of the features...")
            for k in keys:
                report[k] -= np.mean(report[k])

        # get relevant metadata for this tod from pickle file
        tod_name = self.get_name()
        pickle_id = self._pickle_data['name'].index(tod_name)

        # get detectors
        live = store.get(self.inputs.get('dets'))['live_final']
        live_dets = list(np.where(live == 1))[0]

        # store each det timeseries in hdf5
        for tes_det in live_dets:
            if self._downsample:
                data = tod.data[tes_det, ::self._downsample]
            else:
                data = tod.data[tes_det]

            # generate a unique detector id
            det_uid = '%d.%d' % (self.get_id(), tes_det)

            # test if the dataset already exists, if so update it instead
            try:
                dataset = self._group.create_dataset(det_uid, data=data)
            except RuntimeError:
                dataset = self._group[det_uid]
                dataset[:] = data

            # save report to h5 file
            for k in keys:
                dataset.attrs[k] = report[k][tes_det]

            # save label
            dataset.attrs['label'] = int(self._pickle_data['sel'][tes_det, pickle_id])

        self.logger.info("Data saved in %s" % self._output_file)

    def finalize(self):
        self._hf.close()


class PrepareDataLabelNew(Routine):
    def __init__(self, **params):
        """Prepare an HDF5 data set that contains all relevant metedata
        and detector timeseries as a basis for various ML studies
        """
        Routine.__init__(self)
        self.inputs = params.get("inputs", None)
        self._pickle_file = params.get('pickle_file', None)
        self._output_file = params.get('output_file', None)
        self._group_name = params.get('group', None)
        self._remove_mean = params.get('remove_mean', False)
        self._truncate = params.get('truncate', 1000)
        self._downsample = params.get('downsample', 5)

    def initialize(self):
        # load pickle file
        self.logger.info("Loading %s..." % self._pickle_file)
        with open(self._pickle_file, "r") as f:
            self._pickle_data = cPickle.load(f)

        # create output h5 file if it doesn't exist
        if os.path.isfile(self._output_file):
            # file exists
            self.logger.info("Found %s, updating instead" % self._output_file)
            self._hf = h5py.File(self._output_file, 'a')
        else:
            # file doesn't exist
            self.logger.info("Creating %s..." % self._output_file)
            self._hf = h5py.File(self._output_file, 'w')

        try:
            self.logger.info("Creating group %s..." % self._group_name)
            self._group = self._hf.create_group(self._group_name)
        except ValueError:
            self.logger.info("%s exists, update instead" % self._group_name)
            self._group = self._hf[self._group_name]

    def execute(self, store):
        # retrieve tod data
        tod = store.get(self.inputs.get('tod'))
        fft = store.get(self.inputs.get('fft'))['fdata']

        # retrieve the calculated statistics
        report = store.get(self.inputs.get('report'))
        keys = report.keys()

        # remove mean if needed
        if self._remove_mean:
            report = copy.deepcopy(report)
            self.logger.info("Remove means of the features...")
            for k in keys:
                report[k] -= np.mean(report[k])

        # get relevant metadata for this tod from pickle file
        tod_name = self.get_name()
        pickle_id = self._pickle_data['name'].index(tod_name)

        # get list of detectors of interests
        live = store.get(self.inputs.get('dets'))['live_final']
        live_dets = list(np.where(live == 1))[0]

        # treat the median as the common modes
        fdata_cm = np.median(np.abs(fft[live_dets, :self._truncate]), axis=0)

        # store each det timeseries in hdf5
        for tes_det in live_dets:
            tdata = tod.data[tes_det, ::self._downsample][400:400+self._truncate]
            fdata = np.abs(fft[tes_det, :self._truncate])
            fdata -= fdata_cm

            data = np.vstack([tdata, fdata])

            # generate a unique detector id
            det_uid = '%s.%d' % (self.get_name(), tes_det)

            # test if the dataset already exists, if so update it instead
            try:
                dataset = self._group.create_dataset(det_uid, data=data)
            except RuntimeError:
                dataset = self._group[det_uid]
                dataset[:] = data

            # save report to h5 file
            for k in keys:
                dataset.attrs[k] = report[k][tes_det]

            # save label
            dataset.attrs['label'] = int(self._pickle_data['sel'][tes_det, pickle_id])

        self.logger.info("Data saved in %s" % self._output_file)

    def finalize(self):
        self._hf.close()
