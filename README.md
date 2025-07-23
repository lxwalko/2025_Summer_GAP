# LVC MathPhys Research Code
### An attempt at organization

Much of the code is disorganized and undocumented, this GitHub repo is the beginning of an attempt to rectify that.  

### Getting this code onto your machine

You could copy and paste each file into your own text editor, and then save in a location of your choosing. *Don't do this*  
  
Instead, first choose a directory (folder) that you want to save GAP related files in. Remember where this is, and the
filepath to it!  
  
Now, you need to install git. To do so, follow the instructions found [here](https://github.com/git-guides/install-git).  
  
You're now ready to download this code to your machine.
1. Open the terminal (This is called PowerShell or cmd.exe on Windows)
2. Change directories to the directory you created for your GAP files using the command  
`cd file_path_to_your_directory`
3. Run the command `git pull https://github.com/lxwalko/2025_Summer_GAP/tree/master` 

### How to use

*Assuming you have GAP downloaded, installed, and working.*  
    *To do so, visit the [GAP website](https://www.gap-system.org/)*    
    
1. Start GAP from the terminal with the command `gap`
2. Change directory to where this code is saved with the command `ChangeDirectoryCurrent( "filepath" );`
3. Read all the code in this codebase using the command `Read( "load.gap" );`     *Or you can read individual files with* `Read( "filename" );`
  
Documentation (where applicable) is written in the source files. More organized documentation can be found in the `DOC.md` file.
