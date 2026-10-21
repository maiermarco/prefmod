.onUnload <- function(libpath){
  library.dynam.unload(chname = "prefmod", libpath = libpath) # unload prefmod's DLL
}
