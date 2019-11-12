import os,sys
try:
    os.system('python3 untar.py')
except (Exception, e):
    print("Sorry, 'untar.py' has an error.")
    sys.exit(1)

try:
    os.system('python3 grindingdata2.py')
except (Exception, e):
    print("Sorry, 'grindingdata2.py' has an error.")
    sys.exit(1)

try:
    os.system('python3 hist.py')
except (Exception, e):
    print("Sorry, 'hist.py' has an error.")
    sys.exit(1)

try:
    os.system('python3 persistence.py')
except (Exception, e):
    print("Sorry, 'persistence.py' has an error.")
    sys.exit(1)

try:
    os.system('python3 del.py')
except (Exception, e):
    pass
    print("Sorry, 'del.py' has an error.")

