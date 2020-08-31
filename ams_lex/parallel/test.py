import subprocess
for state in range(200, 351, 50):
     for seed in range(10):
             print "N%dP2S%dT32.txt"%(state, seed)
             file = open("N%dP2S%dT32.txt"%(state, seed), "w")
             file.write(subprocess.check_output("~/pure-project/ams_lex/parallel/shortest_ams_lex %d 2 %d 32 0 0"%(state, seed), stderr=subprocess.STDOUT, shell=True))
             file.close()
