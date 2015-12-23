from cluster import *
import optparse
import sys

##########################################
## Options and defaults
##########################################
def getOptions():
    parser = optparse.OptionParser('python *.py [option]"')
    parser.add_option('--sdf',dest='input',help='intput sdf file', default='')
    parser.add_option('--fp',dest='fp',help='fingerprint type: tp,mc,mo (Topological Fingerprints, MACCS Keys, Morgan Fingerprints), default is mc', default='mc')
    parser.add_option('--radius',dest='radius',help=' the radius of the Morgan fingerprint, default is 2',type='int', default=2)   
    parser.add_option('--algorithm',dest='algorithm',help='cluster algorithm :b,m (Butina, Murtagh), default is b', default='b')
    parser.add_option('--cutoff',dest='cutoff',help='distThresh(0-1),elements within this range of each other are considered to be neighbors, needed for Butina cluster algorithm, default is 0.5', type='float', default=0.5)
    parser.add_option('--nclusts',dest='nclusts',help='number of clusters, needed for Murtagh cluster algorithm, default is 1',type='int', default=1)
    parser.add_option('--murtype',dest='Murtype',help='Method for Murtagh:Wards, SLINK, CLINK, UPGMA, needed when Murtagh is set as algorithm, default is Wards', default='Wards')
    parser.add_option('--out',dest='output',help='output sdf file', default='')
    options, args = parser.parse_args()
    
    if options.input=='' or options.output=='':
        parser.print_help()
        print "No input or output is provided"
        sys.exit(1)
    return options

def main():
    options = getOptions()
    fpOpdict = {'tp':'Topological Fingerprints','mc':'MACCS Keys','mo':'Morgan Fingerprints'}
    algOpdict = {'b':'Butina','m':'Murtagh'}
    options.algorithm = algOpdict[options.algorithm]
    print("fingerprint type: %s" % fpOpdict[options.fp])
    if options.fp == 'mo':
        print("radius: %s" % str(options.radius))
    print("cluster algorithm: %s" % options.algorithm)
    if options.algorithm == "Murtagh":
        print("Murtagh method: %s" % options.Murtype)
        print("Murtagh cluster number set: %s" % options.nclusts)
    elif options.algorithm == "Butina":
        print("cutoff(distThresh) : %s" % options.cutoff)
    print 
    print('sdf reading...')
    sdfparse = ChemParse(options.input)
    sdfparse.sdf_reader()
    print('fingerprint calculating...')
    sdfparse.get_fps(options.fp, options.radius)
    print('clustering...')
    fpCluster = Fingerprint_Cluster(sdfparse.fps)
    fpCluster.distance_matrix()
    fpCluster.cluster_dict(options.algorithm, options.cutoff, options.Murtype, options.nclusts)
    print ('done, output to %s' % options.output)
    sdfparse.clusterOutput(options.output, fpCluster.cdict)
    for c in fpCluster.clustdict:
        print("cluster%s: %s" % (str(c),str(fpCluster.clustdict[c])))
    
    
if __name__ == "__main__":
    main()
