import os
import glob

def clean(recipe, ver):
    """Clean up version directory"""
    os.system("rm -v run_{}/errfile_*".format(ver))
    os.system("rm -v run_{}/log_*".format(ver))
    os.system("rm -v run_{}/timefile_*".format(ver))
    os.system("rm -v run_{}/cutparams_*".format(ver))
    os.system("rm -v run_{}/todList_*".format(ver))
    os.system("rm -v run_{}/timefile_*".format(ver))
    # only remove sub-db
    files = glob.glob("run_{}/*.db".format(ver))
    files = [f for f in files if len(os.path.basename(f).split('_'))==6]
    for f in files:
        os.system("rm -v %s" % f.rstrip())

def combine(recipe, ver):
    """Combine db into one"""
    files = glob.glob("run_{}/*.db.*".format(ver))
    for f in files: print f
    files = [f for f in files if os.path.basename(f).split('.')[-1]!='db']
    if len(files) == 0:
        raise Exception("No sub-db files found!")
    filename = '.'.join(files[0].split('.')[:-1])
    mode = "a" if os.path.exists(filename) else "a"
    with open(filename, mode) as ff:
        for f in files:
            with open(f, "r") as tf:
                lines = tf.readlines()
                for l in lines:
                    if l[0] == '#':  # comment
                        continue
                    else:
                        ff.write(l)
