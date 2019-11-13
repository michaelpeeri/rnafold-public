from pyvirtualdisplay.display import Display  # use Xvnc to provide a headless X session (required by ete for plotting)

with Display(backend='xvnc') as disp:  # Plotting requires an X session
    print(disp)
    print("After running this, run 'xhost +' to disable access control...")
    raw_input('Press enter to stop headless X server...')

    
