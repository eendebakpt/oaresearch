import os
import time
from typing import List
from os.path import join

import oapackage

from oapackage import oahelper


class job:

    """ Class representing a job """
    logfile = None
    execute = False

    def __init__(self, cmd, jobtype='generic', shorttag=None, checkfiles=[], checkfilesstart=[], ncores=1, queue=None):
        """ Create an object for job definitions

        Arguments:
            cmd (str): command to be executed
            jobtype (str)
            checkfiles (list, optional): the job is complete when all files are present
            checkfilesstart (list, optional): the job is can be executed if each of the files in the list is present
            ncores (int): number of cores for the job
            queue (str): queue to use in the cluster
        """
        self.cmd = cmd
        self.jobtype = jobtype
        self.verbose = 1
        self.ncores = ncores
        self.queue = queue
        self.shorttag = shorttag
        if isinstance(checkfiles, str):
            self.checkfiles = [checkfiles]
        else:
            self.checkfiles = checkfiles
        if isinstance(checkfilesstart, str):
            self.checkfilesstart = [checkfilesstart]
        else:
            self.checkfilesstart = checkfilesstart

    def writeJobfile(self, outfile):
        f = open(outfile, 'wt')
        f.write('# job: %s' % self.jobtype + os.linesep)
        f.write(self.cmd)
        f.write(os.linesep)
        f.close()

    def __repr__(self):
        c = self.complete()
        if c:
            s = 'job: %s: complete %d' % (self.jobtype, c)
        else:
            s = 'job: %s: complete %d (canrun %d)' % (
                self.jobtype, c, self.canrun())

        return s

    def canrun(self, verbose: int = 0):
        if oahelper.checkFilesOA(self.checkfilesstart, cache=True, verbose=verbose):
            return True
        else:
            return False

    def analyse(jx, verbose=1):
        print(jx)
        if verbose >= 2:
            print('checkfiles:')
            print(jx.checkfiles)
            print('checkfilestart:')
            print(jx.checkfilesstart)
        _ = jx.canrun(verbose=verbose)

    def complete(self, verbose=0):
        if verbose:
            print('job.complete:')
        if oahelper.checkFilesOA(self.checkfiles, verbose=verbose):
            return True
        else:
            return False

    def runjob(self, cache=True, verbose=None):
        """ Run the job """
        goodcalc = 1
        if verbose is None:
            verbose = self.verbose
        if oahelper.checkFilesOA(self.checkfiles, cache=cache):
            if verbose:
                print('job %s: results already calculated' % self.jobtype)
        else:
            if verbose >= 2:
                print('job %s: canrun() %s' % (self.jobtype, self.canrun()))
            if self.canrun():
                if verbose >= 2:
                    print('job %s: running command' % self.jobtype)

                res = oahelper.runcommand(
                    self.cmd, dryrun=(not self.execute), idstr=self.jobtype, logfile=self.logfile)
                if not res == 0:
                    goodcalc = 0

        c = self.complete()
        if self.verbose >= 2:
            print('runjob: c %s, goodcalc %d' % (c, goodcalc))

        if goodcalc == 0:
            c = 0
        return c


def createJobScript(j, index=0, scriptdir=None, queue='q7d', verbose=1, jobtag='quickjob'):
    """ Create a job script file for easy submission

    Arguments:
        j (job object): job to be executed
        index (int): helper variable
        scripdir (str): output directory
        jobtag (str): tag used for output name of the script

    """
    if scriptdir is None:
        basedir = os.path.join(os.getenv('HOME'), 'tmp')
    else:
        basedir = scriptdir
    _ = oapackage.mkdirc(basedir)

    if 'VSC_SCRATCH' in os.environ.keys():
        vsccluster = 1
        if verbose:
            print('we are running on the VSC cluster!')
    else:
        vsccluster = 0

    if queue is None:
        queue = j.queue
    if queue is None:
        queue = 'q24h'

    ncores = j.ncores
    if vsccluster:
        if ncores == 8 or ncores == 6 or ncores == 12 or ncores == 14 or ncores == 18:
            print(
                'warning: jobs on hopper should have number of cores compatible with 20')
            raise Exception('invalid number of cores')

    reservationid = None

    if reservationid is not None:
        resstr = '-l advres=%s' % reservationid
    else:
        resstr = ''

    cmd = j.cmd

    if verbose:
        print(' createJobScript: scriptdir %s' % basedir)

    if j.shorttag is None:
        outfile0 = '%s-%d.sh' % (jobtag, index, )
    else:
        outfile0 = '%s-%d-%s.sh' % (jobtag, index, j.shorttag)
    outfile = os.path.join(basedir, outfile0)

    metafile = join(scriptdir, outfile0.replace('.sh', '-meta.txt'))
    metafile = join(scriptdir, 'meta-information.txt')

    fid = open(outfile, 'w')
    fid.write('#/bin/sh\n\n')
    fid.write('#PBS -l nodes=1:ppn=%d\n' % ncores)
    fid.write('#PBS -q %s' % queue)
    fid.write('# Quick job for VSC %s\n\n' % time.asctime())
    fid.write('# tag: %s' % (jobtag,))

    fid.write('# Job type: %s\n' % j.jobtype)
    fid.write('echo "Job type: %s"\n' % j.jobtype)

    fid.write('DATE=`date`\n')
    fid.write('HOSTNAME=`hostname`\n')
    fid.write('echo "host: ${HOSTNAME}, date: ${DATE}, start" \n')
    fid.write('echo "job: %s" >> %s \n' % (j.jobtype, metafile))
    fid.write('echo "host: ${HOSTNAME}, date: ${DATE}" >> %s \n' % metafile)
    fid.write(
        'echo "   PBS_JOBID $PBS_JOBID, PBS_JOBNAME $PBS_JOBNAME, PBS_QUEUE $PBS_QUEUE" >> %s \n' % metafile)

    fid.write('export OMP_NUM_THREADS=%d;\n' % ncores)

    if vsccluster:
        # for hopper
        if 1:
            fid.write('source /user/antwerpen/201/vsc20149/bin/hopper2017.sh\n')
        else:
            fid.write('module load hopper/2015a\n')
            fid.write('module load CMake/3.1.0-intel-2015a\n')
            fid.write('module load Python/2.7.9-intel-2015a\n')
            fid.write('module load SWIG/3.0.5-intel-2015a-Python-2.7.9\n')
            fid.write('module load matplotlib/1.4.3-intel-2015a-Python-2.7.9')

        fid.write('\n')

    fid.write('\ncd %s\n' % basedir)
    fid.write('%s\n' % (cmd, ))

    fid.write('\necho "quick job %s done"\n' % index)
    fid.write('DATE2=`date`\n')
    fid.write(
        'echo "job %s: date: ${DATE}, done" >> %s \n' % (j.jobtype, metafile))

    fid.write('\n')
    fid.write('\n')

    if vsccluster and (queue is not None) and (queue != 'none'):
        if queue == 'q1h':
            wt = '1:00:00'
        elif queue == 'q24h':
            wt = '24:00:00'
        elif queue == 'q72h':
            wt = '72:00:00'
        else:
            wt = '168:00:00'
        substr = 'qsub %s -q %s -l walltime=%s %s' % (
            resstr, queue, wt, outfile)
    else:
        substr = 'bash  %s' % (outfile)

    fid.write('# submission: %s\n' % (substr))
    fid.close()

    if verbose >= 2:
        print('  writen file %s ' % (outfile))
    if vsccluster:
        if verbose:
            print('  submit with: %s' % substr)

    return outfile, substr


def makeJobList(scriptdir: str, jobs: List, verbose=1, ncores: int = 0, queue=None):
    """ Create a list of jobs

    Args:
        scriptdir (str)
        jobs (list): list of job objects
        verbose (int)
    """
    slist: List = []
    nj = 0
    for idx, j in enumerate(jobs):
        if not j.complete() and j.canrun():
            if ncores > 0:
                j.ncores = ncores
            outfile, substr = createJobScript(
                j, index=idx, verbose=0, queue=queue, scriptdir=scriptdir)
            slist += [substr]
            nj += 1
            if verbose >= 2:
                print('job: %s' % j.shorttag)
            if verbose >= 2 and idx < 5:
                print('  checkfiles: %s' % (j.checkfiles))
            if verbose >= 4:
                j.canrun(verbose=1)
                # j.complete(verbose=1)
                print('  checkfilesstart: %s' % (j.checkfilesstart))
        else:
            if verbose >= 3:
                print('job: %s: complete %d, canrun %d' %
                      (j.shorttag, j.complete(), j.canrun()))
            if verbose >= 3 and idx < 1:
                j.canrun(verbose=1)
                if verbose >= 4:
                    print('  checkfilesstart: %s' % (j.checkfilesstart))
                print('  checkfiles: %s' % (j.checkfiles))

    jfile = join(scriptdir, 'subs.sh')
    if verbose:
        print('created %d/%d jobs in %s' % (nj, len(jobs), jfile))

    fid = open(jfile, 'wt')
    for i, s in enumerate(slist):
        _ = fid.write('%s\n' % s)
    fid.close()
    return jfile
