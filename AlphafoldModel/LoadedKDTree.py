import numpy as np
from scipy.spatial import KDTree
from typing import Callable, Iterable


class LoadedKDTree(KDTree):
   '''
   The `LoadedKDTree` class is an adapter class derived from `scipy.spatial.KDTree` that provides a layer of abstraction to map custom data payloads with spatial information for use with the `scipy.spatial.KDTree` class. It exposes wrapped `scipy.spatial.KDTree` methods to instead return the corresponding data payloads for each coordinate.
   '''
   
   def __init__(self, dataObj: list, coord_getter: Callable[[],Iterable]):
      self.coord_getter = coord_getter
      self.dataObj = dataObj
      self.dataType = type(dataObj[0])
      self.datapoints = [np.array(self.coord_getter(data)) for data in dataObj]
      self.dimensions = len(self.datapoints[0])
      super().__init__(self.datapoints)   # initialize KDTree
      
   
   def query_within_radius(self, query, radius: float, p: float=2, return_sorted: bool=True) -> list:
      '''
      `query: Any` = The query spatial coordinate. Can be an instance of the elements within `dataObj` or a coordinate array.
      
      `radius: float` = The radius within which to search and return neighbouring points.
      
      `p: float` = Which Minkowski p-norm to use. 
      - `1` is the sum-of-absolute-values distance (“Manhattan” distance). 
      - `2` is the usual Euclidean distance. 
      - `infinity` is the maximum-coordinate-difference distance.
      
      `return_sorted: bool` = Return the list of neighbours in order of descending proximity to the query.
      '''
      
      queryCoordinate: list = []
      
      if isinstance(query, self.dataType):
         queryCoordinate = self.coord_getter(query)
      elif isinstance(query, Iterable) and len(query) == self.dimensions:
         queryCoordinate = query
      else:
         error = 'Dimensions of coordinate array does not match existing `dataObj` elements.'
         raise ValueError(error)
      
      neighbours = super().query_ball_point(queryCoordinate, radius, p, return_sorted=return_sorted)
      neighbouringDataObj = [self.dataObj[idx] for idx in neighbours]
      
      return neighbouringDataObj
   
   
   def query_k_nn(self, query, k: int=1, p: float=2, distance_upper_bound: float=float('inf'), inc_distance: bool=False) -> list:
      '''
      Returns a sorted list of k-nearest neighbours from the query.
      
      `query: Any` = The query spatial coordinate. Can be an instance of the elements within `dataObj` or a coordinate array.
      
      `k: int` = The number of nearest neighbours to return.
      
      `p: float` = Which Minkowski p-norm to use. 
      - `1` is the sum-of-absolute-values distance (“Manhattan” distance). 
      - `2` is the usual Euclidean distance. 
      - `infinity` is the maximum-coordinate-difference distance.
      
      `distance_upper_bound: float` = Distance limit from query to search for points.
      '''
      
      queryCoordinate: list = []
      
      if isinstance(query, self.dataType):
         queryCoordinate = self.coord_getter(query)
      elif isinstance(query, Iterable) and len(query) == self.dimensions:
         queryCoordinate = query
      else:
         error = 'Dimensions of coordinate array does not match existing `dataObj` elements.'
         raise ValueError(error)
      
      neighbourDistances, neighbourIdx = super().query(queryCoordinate, k, p=p, distance_upper_bound=distance_upper_bound)
      nearestNeighbourObj = [self.dataObj[idx] for idx in neighbourIdx]
      
      if inc_distance == False:
         return nearestNeighbourObj
      elif inc_distance == True:
         return neighbourDistances, nearestNeighbourObj



if __name__ == "__main__":
   # Sample points in 3D space
   points = np.array([[1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9], [9, 10, 11], [1,30,20], [5,5,5]])

   # Create a KDTree
   tree = KDTree(points)

   # Query all points within a certain radius
   radius = 5.0
   query_point = np.array([4, 5, 6])
   indices = tree.query_ball_point(query_point, radius)

   print("Points within radius", radius, ":", points[indices])
   
   from ModelPDB import alphafold
   
   atoms = alphafold.atoms
   
   def get_coordinates(atom):
      return (atom.x, atom.y, atom.z)
   
   kdtree = LoadedKDTree(atoms,get_coordinates)
   
   neighbours = set([atom.resId for atom in kdtree.query_within_radius(atoms[2000],5)])
   k_nn = kdtree.query_k_nn(atoms[2000],2) 
   print(neighbours)
   print(k_nn)
