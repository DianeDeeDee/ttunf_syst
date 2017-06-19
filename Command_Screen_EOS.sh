screen -S read
pagsh -c /bin/bash
export KRB5CCNAME=FILE:/tmp/diane_krb 
kinit -r 30d
Enter Password
aklog
source ~/.bashrc
source Kinit2.sh
detach screen, 
log into same lxplus machine as before
screen -dr
aklog
source EOS.sh
source RoottufRunLocallySyst.sh

