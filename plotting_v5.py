
import csv
from plotting_functions import *
import os.path

files = os.listdir(src_path)
for file_name in files:
    if file_name[-4:] == '.csv':
        with open(file_name) as csvFile, open(file_name[:-4] + '.txt', "w+") as txtFile:
            reader = csv.reader(csvFile)
            i = 0
            for row in reader:
                if i >= 1:
                    txtFile.write("2019-6-28\t{time}\t{absorbance}\t0\t0\n".format(time=row[0], absorbance=row[1]))
                i += 1
            txtFile.close()
            csvFile.close()

# update ranges for linear fits
u1 = {}
u2 = {}
u3 = {}

### first, change the files names as needed (e.g. 0.3.Sample.csv becomes 0.3.Sample)
### then go down to where all the commands are, "comment" everything (place a # before) except for the allplot() command, then run
### fix the ranges below as necessary, then close all the graphs and now comment allplot, uncomment allplot_linfit and ratevconc
#u1['GOx40.Sample'] = [0.3, 60, 90] #25
#u1['GOx80.Sample'] = [0.5, 60, 90] #1
#u1['1repeat.Sample'] = [1, 45, 70] #20
#u1['1repeat2.Sample'] = [1, 65, 90]
#u1['GOx120.Sample'] = [1, 60, 90]
#u1['GOx160.Sample'] = [2, 60, 90] #19
#u1['3.5repeat.Sample'] = [3.5, 52, 72] #5
#u1['5repeat.Sample'] = [5, 55, 75] #3
#u1['10repeat.Sample'] = [10, 75, 105] #26

#u1['100.Sample'] = [100, 60, 90]
u1['0.3.Sample'] = [0.3, 42, 60]
u1['0.5.Sample'] = [0.5, 41, 59]
#u1['1.Sample'] = [1, 54, 75]
u1['1.Sample'] = [1, 60, 110]
u1['2.Sample'] = [2, 42, 55]
u1['3.5.Sample'] = [3.5, 40, 60]
u1['5.Sample'] = [5, 41, 62]
# u1['10.Sample'] = [10, 40, 80]
IMPD = idata_nooffset(src_path)
fits1 = looplinfit(u1, IMPD)

# fits2 = looplinfit(u2, IMPD)
# fits3 = looplinfit(u3, IMPD)

# for key in needfix:
# L = len(IMPD[0][key][0])
# x = zeros(L)
# for i in range(1,L,2):
# x[i] = -0.5
# IMPD[0][key][0] = IMPD[0][key][0] + x

### these are the commands we use, others are similar:
allplot_linfit(u1, IMPD)
ratevconc(fits1)

### allplot linfit and ratevconc uncommented together
# allplot()

#lb plots may be used later on
lineweaverburk(fits1)
lbcustomfit(fits1)
ratevconc_enzlin(fits3)

