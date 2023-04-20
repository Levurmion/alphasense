from time import time

def function_timer(message):
   
   # function_timer returns the actual wrapper function
   def func_wrapper(func):

      # timed_func takes in the `message` parameter from function_timer
      # `message` is available down to this level via closures
      def timed_func(*args, **kwargs):

         startTime: float = time()
         funcOutput = func(*args, **kwargs)
         endTime: float = time()
         print(message, endTime-startTime, 's')
         return funcOutput
      
      # func_wrapper returns the modified `func`
      return timed_func
   
   # function_timer returns the actual wrapper (func_wrapper)
   #
   # the first decorator is essentially returning the actual decorator
   # with the `message` params made available to it without explicit declaration in `func_wrapper`
   # because of the closure
   return func_wrapper  


# file generator function
def readFile_as_generator(filePath: str):
   for line in open(filePath):
      yield line