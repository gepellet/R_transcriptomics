# adaptative well_name function

Adjust_Well_Name <- function(name_vector,add_0){
  adjusted_name_vector = c()
  for (i in 1:length(name_vector)){#length(name_vector)
    if (add_0 == T){
      index = regexpr(pattern = "_[A-Z][1-9]{1}$",name_vector[i], fixed = F)
      if (length(index) != 0){
        replace = paste(substr(name_vector[i],index,index+1),'0',
                          substr(name_vector[i],index+2,nchar(name_vector[i])),sep='')
        adjusted_name_vector[i] = gsub(pattern = "_[A-Z][1-9]{1}$" ,
                                       replacement = replace ,name_vector[i], fixed = F)
      }
    }else{
      index = regexpr(pattern = "_[A-Z]0[1-9]{1}$",name_vector[i], fixed = F)
      if (length(index) != 0){
        replace = paste(substr(name_vector[i],index,index+1),
                        substr(name_vector[i],index+3,nchar(name_vector[i])),sep='')
        adjusted_name_vector[i] = gsub(pattern = "_[A-Z]0[1-9]{1}$" ,
                                       replacement = replace ,name_vector[i], fixed = F)
      }
    }
  }
  return(adjusted_name_vector)
}