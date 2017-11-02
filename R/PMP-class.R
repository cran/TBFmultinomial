##' Class for PMP objects
##'

setClass(Class="PMP",
         slots= list(posterior="vector",
                     prior="vector",
                     model="vector",
                     method="character",
                     G="vector"))
