# How to use screen:

*Screen* allows you to ssh onto multiple terminal tabs at the same time. This is an advantage because you dont need to ssh into every terminal tab.

  1. ssh to your cluster
  2. in home directory, type *screen* and then Return --> you will now be inside **screen**
  3. C+a then c (press control and a simultaneous, release press c): to create a new window.
  4. C+a then " then arrows: to move to a different tab use arrows up or down and press Return.
  5. C+a then d (Control-a, release, then press d): used to detach from *screen*. After this, you can normally exit from the cluster.
 
Next day:

  6. ssh to your cluster
  7. screen -ls: to look for the screen sessions that you had open last time
  8. screen -r ...: where ... is replaced with the number of the last screen session
  
  
