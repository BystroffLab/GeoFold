# -*- coding: utf-8 -*-
"""
Created on Thu May 28 10:41:31 2015

@author: walcob
"""
#imports
import ransac
import numpy

def getplot(infile):
  first = True
  plot = open(infile,'r')
  for line in plot:
    line = line.split()
    if first:
      data = numpy.array([float(line[0]),float(line[1])])
      x = numpy.array([float(line[0])])
      y = numpy.array([float(line[1])])
      first = False
    else:
      data = numpy.vstack((data,numpy.array([float(line[0]),float(line[1])])))
      x = numpy.vstack((x,numpy.array([float(line[0])])))
      y = numpy.vstack((y,numpy.array([float(line[1])])))
  return data,x,y

def fit(plotfile,outfile):
  data,x,y = getplot(plotfile)
  n = 4
  k = 10000
  t = .5
  d = 0
  debug = False
  input_columns = [0]
  output_columns = [1]
  model = ransac.LinearLeastSquaresModel(input_columns,output_columns,debug=debug)
  ransac_fit, ransac_data = ransac.ransac(data,model,n,k,t,d,debug=True,return_all=True)
  inx = ransac_data['inliers']
  first = True
  for i in inx:
    if first:
      inliers = numpy.array(data[i:i+1,:])
      first = False
    else:
      inliers = numpy.vstack((inliers,numpy.array(data[i:i+1,:])))
  numpy.savetxt(outfile,inliers)


def main():
  data,x,y = getplot("test.plot")
  print(data)
  n = 5
  k = 10000
  t = 1.
  d = 1
  debug = False
  input_columns = [0]
  output_columns = [1]
  model = ransac.LinearLeastSquaresModel(input_columns,output_columns,debug=debug)
  ransac_fit, ransac_data = ransac.ransac(data,model,n,k,t,d,debug=False,return_all=True)
  print(ransac_fit)
  print(ransac_data)
  inx = ransac_data['inliers']
  first = True
  for i in inx:
    if first:
      inliers = numpy.array(data[i:i+1,:])
      first = False
    else:
      inliers = numpy.vstack((inliers,numpy.array(data[i:i+1,:])))
  print(inliers)
  numpy.savetxt("out.plot",inliers)

if __name__ == "__main__":
  main()