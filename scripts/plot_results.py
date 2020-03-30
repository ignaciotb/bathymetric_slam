#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plot
import numpy
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--initial_poses", dest="initial_poses",
                  default="", help="The filename that contains the original poses.")
parser.add_option("--corrupted_poses", dest="corrupted_poses",
                  default="", help="The filename that contains the optimized poses.")
parser.add_option("--optimized_poses", dest="optimized_poses",
                  default="", help="The filename that contains the optimized poses.")
parser.add_option("-e", "--axes_equal", action="store_true", dest="axes_equal",
                  default="", help="Make the plot axes equal.")
parser.add_option("--output_file", dest="outputFile",
                  default="", help="The output file.")

(options, args) = parser.parse_args()

# Read the original and optimized poses files.
poses_original = None
if options.initial_poses != '':
  poses_original = numpy.genfromtxt(options.initial_poses,
                                    usecols = (1, 2, 3))

poses_corrupted = None
if options.corrupted_poses != '':
  poses_corrupted = numpy.genfromtxt(options.corrupted_poses,
                                    usecols = (1, 2, 3))

poses_optimized = None
if options.optimized_poses != '':
  poses_optimized = numpy.genfromtxt(options.optimized_poses,
                                     usecols = (1, 2, 3))

sum_opt = 0.0
sum_corr = 0.0
for i in range(0, len(poses_optimized)):
  sum_opt += numpy.linalg.norm(poses_original[i,:] - poses_optimized[i,:])
  sum_corr += numpy.linalg.norm(poses_original[i,:] - poses_corrupted[i,:])

sum_opt /= len(poses_optimized)
sum_corr /= len(poses_optimized)

print "Diff original and corrupted", sum_corr
print "Diff original and optimized", sum_opt

#  with open(options.outputFile, "a") as text_file:
  #  text_file.write("%s" % sum_corr)
  #  text_file.write(" %s" % sum_opt)
  #  text_file.close()
#
# Plots the results for the specified poses.
figure = plot.figure()

if poses_original is not None:
 axes = plot.subplot(1, 1, 1, projection='3d')

 plot.plot(poses_original[:, 0], poses_original[:, 1], poses_original[:, 2],
           '-', alpha=0.5, color="green", linewidth=2.0)

if poses_corrupted is not None:
 plot.plot(poses_corrupted[:, 0], poses_corrupted[:, 1], poses_corrupted[:, 2],
           '-', alpha=0.5, color="red", linewidth=2.0)

if poses_optimized is not None:
 plot.plot(poses_optimized[:, 0], poses_optimized[:, 1], poses_optimized[:, 2],
           '-', alpha=0.5, color="blue", linewidth=2.0)

plot.title('Trajectories: GT (Green), Corrupted (Red), optimized (Blue)')

plot.show()
