from todloop import Routine
from moby2.scripting import products


class FindRebiasTime(Routine):
    def __init__(self, config_file=None, offset=0, rebias_wait=0,
                 IV_wait=0):
        """A routine to find IV/rebias time. The offset will be saved
        in the store with key 'offset'

        Args:
            config_file (str): manifest_config file (default None)
            offset (int): (default 0)
            rebias_wait (int): (default 0)
            IV_wait (int): (default 0)

        """
        Routine.__init__(self)
        self._config_file = config_file
        self._offset = offset
        self._rebias_wait = rebias_wait
        self._IV_wait = IV_wait

    def execute(self, store):
        # load tod
        obs = self.get_name()
        # Find IV/rebias gap time
        ct = int(obs.split("/")[-1].split(".")[0])
        ctimes = (ct-IV_wait,ct)
        if products._get_instrument() == "actpol":
            from moby2.instruments import actpol as inst
        elif products._get_instrument() == "mbac":
            from moby2.instruments import mbac as inst
        try:
            db = inst.TODDatabase(config_file=self._config_file)
        except:
            db = None
        if (db is None) or (db.db is None):
            self.logger.info("Database not accessible, IV/rebias offsets set to 0")
            offset_IV = 0
            offset_rebias = 0
        else:
            recs = db.select_acqs(suffix='iv', ctime=ctimes)
            if len(recs) > 0:
                offset_IV = (IV_wait - ct + recs[-1].ctime)*400
            else:
                offset_IV = 0
            self.logger.info("IV offset set to %d"%offset_IV)
            ctimes = (ct-self._rebias_wait,ct)
            recs = db.select_acqs(suffix='bc1_step', ctime=ctimes)
            if len(recs) > 0:
                offset_rebias = (self._rebias_wait - ct + recs[-1].ctime)*400
            else:
                offset_rebias = 0
            self.logger.info("Rebias offset set to %d"%offset_rebias)

        offset = max( offset_IV, offset_rebias, self._offset*400 )
        self.logger.info("Total offset set to %d" %offset)

        # save into store
        store.set('offset', offset)
