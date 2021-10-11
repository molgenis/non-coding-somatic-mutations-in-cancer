#!/usr/bin/bash
#https://www.cyberciti.biz/faq/how-to-run-command-or-code-in-parallel-in-bash-shell-under-linux-or-unix/
# Our custom function
cust_func(){
  echo "Do something $1 times..."
  if [ $(($1%2)) -eq 0 ] #https://www.tutorialsandyou.com/bash-shell-scripting/even-odd-14.html
  then
    sleep 1
    echo "Number $1 is even."
  else
    sleep 3
    echo "Number $1 is odd."
  fi
  
}
# For loop 5 times
for i in {1..5}
do
	cust_func $i & # Put a function in the background
done
 
## Put all cust_func in the background and bash 
## would wait until those are completed 
## before displaying all done message
wait 
echo "All done"