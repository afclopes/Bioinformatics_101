# Connecting to cluster:

1. Connect to VPN then:
    1. Connect to cluster by ‘ssh gateway’
    2. Use password
    3. First time you do it: The command-line will tell you: The authenticity of host 'gateway (...)' can't be established. RSA key fingerprint is ….. Are you sure you want to continue connecting (yes/no)? 
    4. Use password
    5. To check if you are connected: >> hostname -f ,reply: >>>gateway.name
    6. To exit the cluster: >> exit

2. Connect without VPN:
    1. Connect to relay with ‘ssh-relay’:  >> ssh username@ssh-relay.address
    2. Use password
    3. Connect to cluster by ‘ssh gateway’
    4. First time you do it: Commandline will tell you: The authenticity of host 'gateway (....)' can't be established. RSA key fingerprint is …… Are you sure you want to continue connecting (yes/no)? 
    5. Use password
    6. To check if you are connected: >> hostname -f ,reply: >>>gateway.name
    7. To exit the cluster: >> exit and again >> exit

Note: you need to ask IT to set this up for you. They need to make an account for you which is different from connecting via VPN.
Note: working on the cluster via connection through VPN are not supposed to be slower than without VPN.

3. Connect without VPN + config file:
    1. Write a config file with the following info:
  Host relay
  Hostname ssh-relay.mrc-mbu.cam.ac.uk
  User username
    2. Save this file in the directory: ~/.ssh
    3. Now, connect to relay by typing one of:
        - >> ssh relay
        - >> ssh ssh-relay.address
        - >> ssh username@ssh-relay.address
        
 4. Connect without VPN + onto lambda
 
 This step is necessary to be able to make plots on the cluster and visualise them on the Mac.
 
 ssh -X username@ssh-relay.address
 
 ssh -X gateway.name.address
 
 All of your work folders should be in there and multiple different packages to use also available without you having to export the path.
 
 ## Different clusters:
 
Some servers on the cluster come with the possibility to connect to the internet, so you can download files from websites (eg. big data files), but it might not have
all the tools you want available by default, so you might have to export the path in this server. Other servers on the cluster might not allow you to connect to websites 
to download files to work on, but these servers might be more up-to-date on the types of tools that you are interested to use. This difference in servers is because they are 
setup differently and one is associated to the public network, while the other is private.
    
  Interesting note: you can go to and from the different servers without having to exit fully one and enter fully the other server. For example, from servername1 (any 
  directory), you can do *ssh servername2* and you will get into the servername2 server. You can then use servername2 to navigate to your desired directory 
  (because regardless of whether you enter servername2 via servername1 or not, you will arrive at the home directory) and there you can download your file from a website.
  Then you can exit servername2 and you will appear again in servername1 in the exact directory where you were before you went to servername2.
  

