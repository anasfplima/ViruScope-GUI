# DECORATOR / WRAPPER FOR THE SCRAPING SCRIPT

# SOURCE: https://stackoverflow.com/questions/21827874/timeout-a-function-windows

import functools
from threading import Thread

#%%

# takes the timeout parameter and returns the decorator
def timeout(sec):
    # takes target function and adds the extra code to it without altering it
    def deco(func):
        @functools.wraps(func)
        
        # replaces func when it's called
        def wrapper(*args, **kwargs):
            # creates list that stores the result of running the function (output or error)
            res = [Exception(f'Function {func.__name__} timeout {sec} seconds exceeded.')]
            
            # helper function that calls func
            def newFunc():
                try:
                    # replaces res[0] with the successful output if it exists
                    res[0] = func(*args, **kwargs)
                except Exception as e:
                    # returns the exception if function fails
                    res[0] = e
            
            # initiates thread that runs newFunc
            t = Thread(target=newFunc)
            # thread ends automatically if main program is exited
            t.daemon = True
            
            try:
                # starts the thread, calls newFunc
                t.start()
                # waits for the thread to finish in {timeout} seconds, otherwise stops waiting
                t.join(sec)
            except Exception as e:
                print('Error starting thread.')
                raise e
                
            # save the result    
            ret = res[0]
            
            # if the result is an Exception, raise it
            if isinstance(ret, BaseException):
                raise ret
            return ret
        return wrapper
    return deco

#%%


