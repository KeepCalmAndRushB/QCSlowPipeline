import os
import pandas as pd
import time
import re


def directory(path):
    """Checks for the existence of reading and writing directory"""
    if not os.path.exists(path):
        print('The directory {}'.format(path) + ' does NOT exist.')
        print('Press enter to exit\n')
        input()
        quit()


def mqdone(df, path, start_delay=15, check_interval=300):
    """Checks if MQ is finished"""
    isrunn = 'Finish writing tables  ' in df['job'].values  # todo: remove spaces at the end
    while not isrunn:
        print('MQ Analysis still not finished. Retrying in 5 minutes.\n')
        time.sleep(check_interval)
        df = pd.read_table(path, sep='\t')
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('/', '').str.replace('(', '').str.replace(')', '')
        isrunn = 'Finish writing tables  ' in df['job'].values
    if isrunn:
        print('MQ Analysis is finished. So, starting in {} seconds...\n'.format(start_delay))
        time.sleep(start_delay)
        return isrunn


def peakproperties(df):
    """Checks if MQ anaysis was done with Peak Properties ON"""
    ispeak = list(df['value'][df['parameter'] == 'Calculate peak properties'] == 'True')[0]
    if not ispeak:
        print('Anaysis done without Peak Properties ON!!!\nRe-run Maxquant with correct settings.\nPress enter to exit')
        input()
        quit()
    else:
        return ispeak


def experiments(df):
    """Checks if experiment names were defined"""
    isexpe = True if 'experiment' in list(df.columns) else False
    if not isexpe:
        print('No Experiment names defined!!!\nRe-run Maxquant with correct settings.\nPress enter to exit\n')
        input()
        quit()
    else:
        return isexpe


def nrawsandexper(lraw, lexp, verbose=False):
    """Checks if number of experiments equals the number of raw files"""
    n_raws = len(lraw)
    n_expe = len(lexp)
    if verbose:
        print('We have {} raw files and {} experiments.'.format(n_raws, n_expe))
    if n_raws != n_expe:
        print('n Experiments different from n Rawfiles!!!\nRe-run Maxquant with correct settings.\nPress enter to exit\n')
        input()
        quit()
    else:
        if verbose:
            print('All Ok!\n')


def uniqueexper(dictionary, n, verbose=False):
    """Checks if each experiments is associated to a single raw file"""
    if len([k for k, v in dictionary.items() if list(dictionary.values()).count(v) == 1]) == n:
        if verbose:
            print("{" + "\n".join("{}: {}".format(k, v) for k, v in dictionary.items()) + "}")
            print('\nEach experiments is associated to only one raw file. All OK until now.\n')
    else:
        print('\nYou have at least one experiments associated to more than one raw file!!!\nRe-run Maxquant with correct settings.\nPress enter to exit')
        input()
        quit()


def lfqison(df):   # todo: fix this horror
    """Checks the mqpar.xml file to see if Lfq is ON"""
    regexp = re.compile('(.*?)(<lfqMode>)([0-9])(</lfqMode>)(.*?)') 
    lfqstring = ''
    mqpar = open(df)
    for line in mqpar:
        match = regexp.match(line)
        if match:
            lfqstring = ('%s' % (match.group(3)))
            islfq = bool(int('%s' % (match.group(3))))
            mqpar.close()
            return islfq
    if not bool(lfqstring):
        mqpar.close()
        print('Something went wrong. I didn\'t find <lfqMode>... \nPress enter to exit\n')
        input()
        quit()


def writingdirectory(path, label):
    """Checks for the existence of final writing directory
    and it usually creates that"""
    if not os.path.exists(os.path.join(path, label)):
        os.makedirs(os.path.join(path, label))
    return path + label + '/' + label + '_'
