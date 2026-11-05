# --- Added on Nov 5, 2026
# --- Functions to get statistics from sephora´s phenoparams

phenoDatesToDoY <- function(start = 1, end = 12, phenoDates, 
                            totalDoYs = c(0, cumsum(c(31,28,31,30,31,30,31,31,30,31,30,31))) ){
  
  nameMonth <- c("Jan.", "Feb.", "March", "April", "May", "June", "July",
                 "Au.", "Sep.", "Oc.", "No.", "De.")
  
  MoY <- numeric(length(phenoDates))
  DoM <- numeric(length(phenoDates))
  period <- start:end
  for( i in 1:length(MoY) ){
    MoY[i] <- period[ sum( totalDoYs - phenoDates[i] < 0 ) ]
    
    DoM[i] <- phenoDates[i] - totalDoYs[MoY[i]] 
  }
  
  as.data.frame(cbind(month=nameMonth[MoY], day=DoM))
}

redondea <- function(x){
  ifelse( x %% 1 < 0.5, floor(x), ceiling(x))
}

sephora_basicStats <- function(data){
  x <- vector("list",6)
  
  x[[1]] <- c(get_basicStats(x=data[,1])$median,get_basicStats(x=data[,1])$mad)
  x[[2]] <- c(get_basicStats(x=data[,2])$median,get_basicStats(x=data[,2])$mad)
  x[[3]] <- c(get_basicStats(x=data[,3])$median,get_basicStats(x=data[,3])$mad)
  x[[4]] <- c(get_basicStats(x=data[,4])$median,get_basicStats(x=data[,4])$mad)
  x[[5]] <- c(get_basicStats(x=data[,5])$median,get_basicStats(x=data[,5])$mad)
  x[[6]] <- c(get_basicStats(x=data[,6])$median,get_basicStats(x=data[,6])$mad)
  x
}

get_basicStats <- function(x){
  temp <- getTheStats(x=x)
  
  list(median=temp$dates[3],mad=temp$dates[8], 
       median_text=temp$text[3],mad_text=temp$dates[8])
}

getTheStats <- function(x){
  x <- x[!is.na(x)]
  
  TEMP <- summary(x, na.rm = TRUE)
  NAMES <- names(TEMP)
  SD <- sd(x, na.rm = TRUE)
  MAD <- mad(x, na.rm = TRUE)
  
  df_text <- c()
  for(i in 1:length(TEMP)){
    df_text[i] <- paste("<b> ", NAMES[i], ": ", "</b>", "<font color=\"red\">", 
                        # round(TEMP[i],3), "</font>")
                        ifelse(TEMP[i]%%1 < 0.5, floor(TEMP[i]), ceiling(TEMP[i])), 
                        "</font>")
  }  
  
  df_text <- c(df_text, paste("<b> Std. Dev.:</b>", "<font color=\"red\">",
                              round(SD,3), "</font>"))
  # ifelse(SD%%1 < 0.5, floor(SD), ceiling(TEMP[i])), , "</font>"))
  
  df_text <- c(df_text, paste("<b> M.A.D.:</b>", "<font color=\"red\">", 
                              MAD, "</font>"))
  
  list(text=df_text, dates=c(round(TEMP), Std=round(SD,3), mad=MAD), size=length(x))
}
