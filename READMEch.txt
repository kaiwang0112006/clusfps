一. molfps所需第三方程序安装：

(一). RDKit安装：
1. 复制RDKit_2015_03_1到d:/RDKit_2015_03_1
2. 添加环境变量：RBASE：d:/RDKit_2015_03_1
3. 添加环境变量：PYTHONPATH：%RBASE%
4. 添加环境变量：PATH：%RDBASE%/lib

二. 聚类原理

molfps是对RDKit的包装，以sdf为输入，包括计算fingerprint和聚类两部分，fingerprint包括Topological Fingerprints, MACCS Keys, Morgan Fingerprints三类，Morgan Fingerprints与ECFP和FCFP类似，因此这类fingerprint需要制定中心半径（radius），聚类算法包括Butina, Murtagh两种，Butina适合大量化合物的聚类，算法效率高，需要预先指定类最多成员数(cutoff/distThresh)，Murtagh是层次聚类，需要预先指定分类数

具体例子见example.txt
