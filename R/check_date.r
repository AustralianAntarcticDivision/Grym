
#' @description
#' \lifecycle{experimental}
#' In order to check the dates when using 365 increments in a Grym simulation this function helps translate the dates into increments or the other way around. The function will warn you if either of the dates fall in a leap year and the increment number or date may be out by 1. If this occurs, recommend changing the years to get correct increment numbers.
#'
#' @title Checking increment dates
#'
#' @param inc_1 The start date of the projection year in the ISO8601 format "YYYY-MM-DD"
#' @param inc_2 Either an increment number or a date in the ISO8601 format "YYYY-MM-DD"
#'
#' @return When inc_2 is an increment number `check_dates` returns the date that increment refers to. If given a date `check_dates` returns the increment number in the year from the start date.
#' @export
#'
#' @examples
#' check_dates(inc_1="2018-12-01", inc_2="2019-01-21")
#' check_dates(inc_1="2018-12-01", inc_2=94)
#'

check_dates<-function(inc_1="2019-01-01", inc_2=1){
  if(missing(inc_1)){stop("inc_1 is missing and must be a date in the format of YYYY-MM-DD")}
  leapyears<-seq(from=1900, to=2200, by= 4)

  if(!missing(inc_1) && !inherits(inc_1,"Date")) {
    inc_1<-as.Date(inc_1)
  }

  if(!missing(inc_1) && inherits(inc_1,"Date")) {
    #leapyears<-seq(from=1900, to=2048, by= 4)
    year<-as.numeric(format(inc_1,'%Y'))
    if(year %in% leapyears) message("Warning: inc_1 is in a leap year so inc_2 date may be a day out if using 365 increments in the Grym.")
  }

  if(inherits(inc_2,"numeric")){
    inc2date<-inc_1 + (inc_2-1) #minus one to account for the 1st increment.
    year2<-as.numeric(format(inc2date,'%Y'))
    if(year2 %in% leapyears){
      message("Warning: inc_2 is in a leap year so it may be a day out if using 365 increments in the Grym.")}
    msg<-paste("Increment", inc_2, "from the start date", inc_1, "is:", inc2date, sep=" ")
    message(msg)}

  if(inherits(inc_2,"Date")||inherits(inc_2,"character")){
    if(inherits(inc_2,"character")){inc_2<-as.Date(inc_2)}
    year2<-as.numeric(format(inc_2,'%Y'))
    if(year2 %in% leapyears){
      message("Warning: inc_2 is in a leap year so it may be a day out if using 365 increments in the Grym.")}
    increment2<-as.numeric(inc_2-inc_1)+1
    msg<-paste("Date:", inc_2, "is increment", increment2, "when the first increment is:",inc_1, sep=" ")
    message(msg)}
}
