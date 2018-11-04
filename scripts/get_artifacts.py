#!/bin/python
#
# Get appveyour artifacts

from coreapi import Client
import urllib.request
import os
import platform

#%% Get information on latest build

if platform.system()=='Windows':
    targetdir = r'c:\\svn\\oapackage\\dist'
else:
    targetdir = '/home/eendebakpt/misc/oa/oacode/dist/'

username = 'eendebakpt'
project_name = 'oapackage-4lws8'
baseurl = 'https://ci.appveyor.com/api'

client = Client()
document = client.get(baseurl+'/projects/%s/%s/history?recordsNumber=10' % (username, project_name))

latest_build = document['builds'][0]
build = client.get(baseurl+'/projects/%s/%s/builds/%d' % (username, project_name, latest_build['buildId']))


#%% Get artifacts

for jobtag, job in enumerate(build['build']['jobs']):
    print('job %s: %s' % (job['jobId'], job['name']))

    url = '/buildjobs/{0}/artifacts'.format(job['jobId'])
    artifacts = client.get(baseurl + url)
    print('found %d artifacts' % len(artifacts))
    for artifact in artifacts:
        filename0 = artifact['fileName'].split('/')[-1]
        url = '/buildjobs/{0}/artifacts/{1}'.format(job['jobId'], artifact['fileName'])
        print('   ' + filename0)

        urllib.request.urlretrieve(baseurl + url, os.path.join(targetdir, filename0))
