# How to transfer a file from your computer to the cluster:

Transfer your file from the computer to the relay, then from the relay to the cluster
1. From your computer do 
      >> scp oocytes_Bulk_config.sh username@ssh-relay.address:.
      
      -> it will ask you for your password
      
      -> it will show you when the full file has been transferred
      
      - scp is a type of copy to a remove system
      - : is used to separate the name of your computer from the name of the folder you want to put your files into the cluster
      - . means to copy to the home directory. You could have instead specified the file path to the folder that you wanted
      - remember that there are quicker ways to transfer more files at a time using for example ESC_* to get all files starting with ESC_ or *.sh to
      get all files ending with .sh
      
2. Now you must enter the ssh-relay, and then continue in transfering the files from this relay to the cluster: 
      >> scp *.sh servername:.
      
      -> it will ask for your password
      
      -> it will show you when the full file has been transferred
      
      - you can only see if the files are really there, once you enter the cluster.
      
3. Once in the cluster, then move your files to the folder needed >> mv *.sh foldername/

4. Good to remember: it is never a problem to have many terminals connected to the cluster at the same time.

# How to transfer a file from the cluster to your computer:

Transfer your file from the cluster to the relay, then from the relay to the computer
1. From the cluster do
      >> scp *.html username@ssh-relay.address:.
      
      -> it will ask you for your password
      
      -> it will show you when the full file has been transferred
      
      - what you are saying here is: do a secure copy of all files ending in .html in the directory that you are in and this is going to the address of 
      username@ssh-relay.address, and to the home folder (.) of that address
      
2. Now you must transfer from the relay to the computer. But, this works the other way around!! You must go to a terminal on your computer, and then request for 
those files, like this:
      >> scp username@ssh-relay.address:*.html .
      - Here, what you are saying is to: do a secure copy from the address username@ssh-relay.address, and the copy will be of all files ending in .html and 
      going to the home folder of where you are (you computer). You can add instead of . the address in your computer (folder path) to where you want to store 
      the files. Or better, you can move to your desired folder, then do the transfer, and using . it will know that you want to transfer the files to the current 
      directory.
      
      
      
