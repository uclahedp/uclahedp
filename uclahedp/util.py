import os
import psutil
import time
import numpy as np

def mem():
    process = psutil.Process(os.getpid())
    print(str(round(process.memory_info().rss*1e-6,2)) + ' MB') 
    

def timeTest(t0=None):
    if t0 is None:
        return time.time()
    else:
        deltaT = time.time() - t0
        print('Execution time: '  + str( round(deltaT,3 )  ) + ' s' )
        return deltaT
    
def timeFormat(t):
    sec = t % 60.0
    minutes = round(t/60)
    return "{:02.0f}:{:02.0f}".format(minutes, sec)


     

class timeRemaining():
    """The timeRemaining class provides progress updates and a prediction of time
    remaining in a loop based on previous performance. To use, intialize the
    object outside the loop with the number of iterations as the argument:
        tr = util.timeRemaining(nloops)
    then update it each iteration with the current iteration number.
        tr.updateTimeRemaining(i)
    The object will automatically keep track of the time between each call,
    and use that to print out the estimated time until completion. 
    Printouts are printed every reportevery iterations. 
    """
    def __init__(self, nsteps, reportevery=100):
        self.last_time = time.time()
        self.time_per_step = []
        self.nsteps = nsteps
        self.reportevery = reportevery

    def updateTimeRemaining(self, i):
        this_time = time.time()
        self.time_per_step.append(this_time-self.last_time)
        self.last_time = this_time
        tremain = np.mean(self.time_per_step)*(self.nsteps - i)
        
        if i > 0 and i % self.reportevery == 0:
            print(str(i) + '/' + str(self.nsteps) + ' done ' + 
                              timeFormat(tremain) + ' remaining' )
        
        
    def avgExecutionTime(self):
        return np.mean(self.time_per_step)



if __name__ == "__main__":
    print(timeFormat(92.123))