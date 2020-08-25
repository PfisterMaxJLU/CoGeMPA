import subprocess
import sys

def auto_install_requirements():

    requirements_list = ["requests == 2.24.0", "XlsxWriter == 1.2.7", "matplotlib == 3.3.0", "numpy == 1.19.1", "appJar == 0.94.0"]

    for elem in requirements_list:

        subprocess.check_call([sys.executable, "-m", "pip", "install", elem])

try:
    print("If the problem is with kofamscan and its requirements, this script wont help. Here is a tutorial for that: https://taylorreiter.github.io/2019-05-11-kofamscan/")
    auto_install_requirements()
    
except:
    print("If this script does not work, add pip and python to PATH. Here is a tutorial: https://www.youtube.com/watch?v=UTUlp6L2zkw")





