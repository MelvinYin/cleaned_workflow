For ec2 instance, remember to:
sudo yum install cpan
(also the sudo install list in docker readme), and probably perl. 
sudo cpan::Parser
sudo cpan::XML (or something like that)
sudo cpan::Simple
set and assign security group. 

These are necessary to build meme. To check if meme is built correctly, make sure there is a folder called libexec in meme, with the files (at least) mast_xml_to_txt and mast_xml_to_html. For mast in static, scp from meme/bin. Note that the mast binary must be built from ec2, if it comes from outside it won't be able to find (for some reason) mast_xml_to_txt. For that, it is probably due to pathing, somewhere in mast.c or other scripts, find using the error message, it's a script with 3 directories only. 

When scp-ing the UI folder in aws, remember to copy-paste the old meme folder, and the mast binary from the previous UI folder, into static.

The static folder is necessary for bokeh to work properly. main.py need to be named as such for the bokeh serve UI to work. 

Commands:

(at dir with key) ssh -i "AWS_Keypair_15032019.pem" ec2-user@3.18.150.48
scp -i "AWS_Keypair_15032019.pem" * ec2-user@3.18.150.48:UI
nohup bokeh serve UI --allow-websocket-origin=ec2-3-18-150-48.us-east-2.compute.amazonaws.com:5001 --port 5001
followed by ctrl+z and then bg, to run in background. 

Remove nohup to see output.

To kill, 
ps aux | grep -i bokeh
then kill 15936, or whatever pid is listed from above. 
