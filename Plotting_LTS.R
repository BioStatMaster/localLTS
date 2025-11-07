
########## packages #########
#' Attempt to load an optional dependency without installing it.
#'
#' @param package Package name as a string.
#' @param attach Whether to call `library()` when the package is available.
load_optional_package <- function(package, attach = TRUE) {
  if (requireNamespace(package, quietly = TRUE)) {
    if (attach) {
      library(package, character.only = TRUE)
    }
    TRUE
  } else {
    warning(sprintf("Package '%s' is not available. Related functionality may be skipped.", package))
    FALSE
  }
}

load_optional_package("ggplot2")
load_optional_package("gridExtra")

################# Metric definitions ###############
evaluation_metric_names <- c(
  "precisition", "recall", "precisition_edge", "recall_edge",
  "true_positives", "true_negatives", "false_positives", "false_negatives",
  "true_positives_edge", "true_negatives_edge", "false_positives_edge", "false_negatives_edge",
  "hamming_dist", "hamming_dist_edge", "hamming_dist_ordering",
  "true_graph_total_edges", "estimated_graph_total_edges"
)

metric_key_from_index <- function(r_index) {
  switch(as.character(r_index),
         "1" = "DAG",
         "2" = "MEC",
         "3" = "Oracle",
         stop("r_index must be 1 (DAG), 2 (MEC), or 3 (Oracle)"))
}

extract_metric_matrix <- function(records, metric_key) {
  if (length(records) == 0) {
    return(matrix(numeric(0), nrow = length(evaluation_metric_names), ncol = 0))
  }
  do.call(cbind, lapply(records, function(rec) {
    values <- rec[[metric_key]]
    if (is.null(values)) {
      rep(NA_real_, length(evaluation_metric_names))
    } else {
      as.numeric(values)
    }
  }))
}

############## Result directory helpers #############
default_result_root <- file.path("D:/2022LTS_GSEM", "Result")
if (!exists("result_root", inherits = FALSE)) {
  if (dir.exists(default_result_root)) {
    result_root <- default_result_root
  } else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && !is.null(ctx$path) && nzchar(ctx$path)) {
      result_root <- dirname(ctx$path)
    } else {
      result_root <- getwd()
    }
  } else {
    result_root <- getwd()
  }
}

################# simulation_data_fun ###############
simulation_data_fun = function(N, P, d = 1, method, r_index = 1, num_out = 0, h_ratio = 1, thresh = 2, result_dir = getwd() ){
  result_precision = result_recall = result_hamming = result_implementations = NULL
  result_precision_byP = result_recall_byP = result_implementations_byP = NULL
  graph_N = rep(0, length(N))
  MRS_N = MRS_N_SE = matrix(0, nrow = length(N), ncol = length(evaluation_metric_names))
  MRS_N_implementation = matrix(0, nrow = length(N), ncol = 1)
  tryCatch({
    for(p_index in 1:length(P)){
      for(n_index in 1:length(N)){
        p_real = P[p_index]; n_real = N[n_index]
        if(method =="LTS") DAGresult_filename = paste0("LTS_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_h",h_ratio * 100,"_thresh",thresh,"_result.Rdata")
        if(method =="USB") DAGresult_filename = paste0("USB_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        if(method =="HGSM") DAGresult_filename = paste0("HGSM_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        if(method =="HLSM") DAGresult_filename = paste0("HLSM_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        # HSLM_p
        if(method =="GDS") DAGresult_filename = paste0("GDS_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        if(method =="TD") DAGresult_filename = paste0("TD_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        full_path = file.path(result_dir, DAGresult_filename)
        MRS_out = data.loading.fun(full_path, r_index)
        MRS_N[n_index, ] = MRS_out$mean_output
        MRS_N_SE[n_index, ] = MRS_out$se_output
        MRS_N_implementation[n_index, ] = MRS_out$number_of_implementation
        graph_N[n_index] = data.loading.fun.ham(full_path, r_index = 1)
      }
      result_precision_byP =  data.frame(x= N,y = MRS_N[,1], SE = MRS_N_SE[,1],
                                         Algorithm = rep( paste0(method),length(N) ) )
      result_precision = rbind(result_precision, result_precision_byP)
      result_recall_byP =  data.frame(x= N,y = MRS_N[,2], SE = MRS_N_SE[,2],
                                      Algorithm = rep( paste0(method),length(N) ) )
      result_recall = rbind(result_recall, result_recall_byP)
      result_hamming_byP =  data.frame(x= N,y = MRS_N[,13], SE = MRS_N_SE[,13],
                                       Algorithm = rep( paste0(method),length(N) ) )
      result_hamming = rbind(result_hamming, result_hamming_byP)

      result_implementations_byP =  data.frame(x= N,y = MRS_N_implementation,
                                               Algorithm = rep( paste0(method),length(N) ) )
      result_implementations = rbind( result_implementations, result_implementations_byP )
    }
  }, error = function(e) {
    return(paste0("the error code: '", e, "'"))
  })
  return( list(result_precision = result_precision, result_recall = result_recall, result_implementations = result_implementations, result_hamming = result_hamming, graph_precision = graph_N ) )
}

################# data.loading.fun ###############
data.loading.fun = function(DAGresult_filename, r_index = 1){

  if (!file.exists(DAGresult_filename)) {
    stop(sprintf("Result file '%s' not found", DAGresult_filename))
  }
  load(DAGresult_filename)

  if (!exists("evaluation_result")) {
    stop("evaluation_result object not found in Rdata file")
  }

  metric_key <- metric_key_from_index(r_index)
  metric_matrix <- extract_metric_matrix(evaluation_result, metric_key)

  if (!is.null(metric_matrix) && length(metric_matrix) > 0) {
    mean_output  = rowMeans(metric_matrix, na.rm = TRUE)
    se_output  = apply(metric_matrix, 1, sd, na.rm = TRUE)/sqrt(ncol(metric_matrix))
    number_of_implementation = ncol(metric_matrix)
  } else {
    mean_output = rep(NA_real_, length(evaluation_metric_names))
    se_output = rep(NA_real_, length(evaluation_metric_names))
    number_of_implementation = 0
  }

  return( list(mean_output = mean_output, se_output = se_output, number_of_implementation = number_of_implementation) )
}

################# data.loading.fun ###############
data.loading.fun.ham = function(DAGresult_filename, r_index = 1){

  if (!file.exists(DAGresult_filename)) {
    stop(sprintf("Result file '%s' not found", DAGresult_filename))
  }
  load(DAGresult_filename)

  if (!exists("evaluation_result")) {
    stop("evaluation_result object not found in Rdata file")
  }

  metric_key <- metric_key_from_index(r_index)
  metric_matrix <- extract_metric_matrix(evaluation_result, metric_key)
  if (is.null(metric_matrix) || ncol(metric_matrix) == 0) {
    return(0)
  }
  hamming_values <- metric_matrix[13, , drop = TRUE]
  sum(hamming_values == 0, na.rm = TRUE)
}

############## Complete Result #################
# Select the scenario folder containing the `.Rdata` files saved by the
# simulation script. Adjust `current_model` as needed (e.g. "OutlierAll",
# "OutlierAll_fraction100").
current_model = "Outlier1"
result_directory = file.path(result_root, current_model)
if (!dir.exists(result_directory)) {
  stop(sprintf("Result directory '%s' not found.", result_directory))
}

N = seq(50, 500, by = 50)
#N = c(200, 500, 1000)
d = 1
P = 5
j = 1
r_index = 4;
num_out = 1
### j: 1 DAG, 2 MEC, 3: Oracles ###
### r_index: 1 precision, 2 recall, 3 implementation, 4 hamming ####
################### plotting ############
result1 = result2 = result3 = result4 = result5 = result6 = result7 = NULL
  #result1 = simulation_data_fun(N, P, d, "HGSM", j)[[r_index]]
  result1 = simulation_data_fun(N, P, d = 1, "LTS", j, num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = result_directory)[[r_index]]
  result2 = simulation_data_fun(N, P, d = 1, "LTS", j, num_out = num_out, h_ratio = 0.7, thresh = 2, result_dir = result_directory)[[r_index]]
  result3 = simulation_data_fun(N, P, d = 1, "LTS", j, num_out = num_out, h_ratio = 0.9, thresh = 2, result_dir = result_directory)[[r_index]]

  result4 = simulation_data_fun(N, P, d, "USB", j, num_out = num_out, result_dir = result_directory)[[r_index]]
  result5 = simulation_data_fun(N, P, d, "GDS", j, num_out = num_out, result_dir = result_directory)[[r_index]]
  result6 = simulation_data_fun(N, P, d, "HLSM", j, num_out = num_out, result_dir = result_directory)[[r_index]]
  result7 = simulation_data_fun(N, P, d, "TD", j, num_out = num_out, result_dir = result_directory)[[r_index]]
  
  result1[,4] = rep("LTS0.5",length(N))
  result2[,4] = rep("LTS0.7",length(N))
  result3[,4] = rep("LTS0.9",length(N))
  
  result = as.data.frame( rbind( result1, result2,result3, result4, result5, result6, result7) )

  result[,4] <- factor(result[,4], levels = c("LTS0.5","LTS0.7", "LTS0.9", "USB", "HLSM") )
  
linesize= 1.0
#id = 1:5
#id = c(1,2,3,4,5,6)
id = c(1,2,3,4,5,6,7)
Color = c("black", "red", "blue", "purple4", "green4", "tomato4", "orange")[id]
LineType = c("solid", "solid", "dashed", "dashed",  "longdash" , "longdash",  "dotted")[id]
LineShape = c(16, 17, 18, 6, 7, 8, 12)[id]

### png 450/300 ###
q = ggplot( data = result, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
q = q + geom_errorbar(data = result, aes(ymin = y-SE ,  ymax = y+SE), width= 10)
q = q + geom_line( aes(linetype=Algorithm), size = linesize ) +
  scale_color_manual(values=Color  ) +
  scale_linetype_manual(values= LineType) +
  scale_shape_manual(values = LineShape)
q = q + geom_point( size = 2.5 )
if(r_index < 3){
  q = q + xlab("Sample Size") + ylab("")  + ylim(0.0,1)
}else{
  q = q + xlab("Sample Size") + ylab("") + ylim(0.0,50)
}
q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "right" )
print(q)

#################################################
result = cbind( result1[,2], result2[,2],result3[,2], result4[,2], result5[,2])
if (load_optional_package("xtable")) {
  print(xtable(result))
} else {
  warning("Package 'xtable' not available; skipping LaTeX table output.")
}

model = c('OutlierAll', 'Outlier1')

N = seq(50, 500, by = 50)
###### plots for comparison of algorihthms ### ####
for(model in c("OutlierAll")){
  model_dir = file.path(result_root, model)
  if (!dir.exists(model_dir)) {
    warning(sprintf("Result directory '%s' not found; skipping.", model_dir))
    next
  }
  for(P in c(20)  ){
    for( d in c(1) ){
      for(num_out in c(1)){
        for(r_index in c(1,2,4) ){
          if(r_index ==1){
            plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"P",P,"D",d,"B",num_out,"PrecisionDAG.png")
          }else if(r_index ==2){
            plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"P",P,"D",d,"B",num_out,"RecallDAG.png")
          }else if(r_index ==4){
            plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"P",P,"D",d,"B",num_out,"HamDAG.png")
          }
          png(plot_filename, width = 450, height = 300)
  
          ################### plotting ############
          result1 = result2 = result3 = result4 = result5 = result6 = result7 = NULL
            result1 = simulation_data_fun(N, P, d, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[r_index]]
            result2 = simulation_data_fun(N, P, d, "HLSM", num_out = num_out, result_dir = model_dir)[[r_index]]
            result3 = simulation_data_fun(N, P, d, "HGSM", num_out = num_out, result_dir = model_dir)[[r_index]]
            result4 = simulation_data_fun(N, P, d, "TD", num_out = num_out, result_dir = model_dir)[[r_index]]
            result5 = simulation_data_fun(N, P, d, "USB", num_out = num_out, result_dir = model_dir)[[r_index]]
            result6 = simulation_data_fun(N, P, d, "GDS", num_out = num_out, result_dir = model_dir)[[r_index]]
          result1[,4] = rep("RGSM", length(N))
          result5[,4] = rep("US", length(N))
          result = as.data.frame( rbind( result1, result2, result3, result4, result5, result6) )
          
          result[,4] <- factor(result[,4], levels = c("RGSM","HLSM", "HGSM", "TD","US", "GDS"))
  
          linesize= 1.0
          Color = c("black", "red", "blue", "purple4", "green4", "tomato4", "orange")[id]
          LineType = c("solid", "solid", "dashed", "dashed",  "longdash" , "longdash",  "dotted")
          LineShape = c(16, 17, 18, 6, 7, 8, 12)
  
          ### png 450/300 ###
          q = ggplot( data = result, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
          q = q + geom_errorbar(data = result, aes(ymin = y-SE ,  ymax = y+SE), width= 10)
          q = q + geom_line( aes(linetype=Algorithm), size = linesize ) +  
            scale_color_manual(values=Color) +
            scale_linetype_manual(values= LineType) +
            scale_shape_manual(values = LineShape)
          q = q + geom_point( size = 2.5 ) 
          if(r_index < 3){
            q = q + xlab("Sample Size") + ylab("")  + ylim(0.0,1)
          }else{
            if(P <150){
              q = q + xlab("Sample Size") + ylab("") # + ylim(0.0,1000)
            }else{
              q = q + xlab("Sample Size") + ylab("") + ylim(0,100)
            }
          }
          q = q + theme_minimal(base_size = 17, base_family = "",
                                base_line_size = 11/22, base_rect_size = 10/100)
          q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
          print(q) 
          
          dev.off()
          # Sys.sleep(1)
}}}}}







############# LTS Performce by P ###########
N = seq(50, 500, by = 50)
N2 = N/2
scale_param = TRUE
scale_param = FALSE
for(model in c("OutlierAll")){
  for(num_out in c(1)){
    model_dir = file.path(result_root, model)
    if (!dir.exists(model_dir)) {
      warning(sprintf("Result directory '%s' not found; skipping.", model_dir))
      next
    }
    if(scale_param){
      plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"D",d,"B",num_out,"C.png")
    }else{
      plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"D",d,"B",num_out,".png")
    }
    png(plot_filename, width = 450, height = 300)
  
    ################### plotting ############
    result1 = result2 = result3 = result4 = result5 = result6 = result7 = NULL
    result001 = simulation_data_fun(N, P = 5, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result002 = simulation_data_fun(N, P = 10, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result003 = simulation_data_fun(N, P = 15, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result004 = simulation_data_fun(N, P = 30, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result005 = simulation_data_fun(N, P = 50, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
  
    if( scale_param ){
      result1 = data.frame( x = N2/log(5) , y = result001, SE = 0, Algorithm = rep("P=5",length(N))  )
      result2 = data.frame( x = N2/log(10), y = result002, SE = 0, Algorithm = rep("P=10",length(N))  )
      result3 = data.frame( x = N2/log(15), y = result003, SE = 0, Algorithm = rep("P=15",length(N))  )
      result4 = data.frame( x = N2/log(20), y = result004, SE = 0, Algorithm = rep("P=20",length(N))  )
      result5 = data.frame( x = N2/log(25), y = result005, SE = 0, Algorithm = rep("P=25",length(N))  )
    }else{
      result1 = data.frame( x = N2, y = result001, SE = 0, Algorithm = rep("P=5",length(N))  )
      result2 = data.frame( x = N2, y = result002, SE = 0, Algorithm = rep("P=10",length(N))  )
      result3 = data.frame( x = N2, y = result003, SE = 0, Algorithm = rep("P=15",length(N))  )
      result4 = data.frame( x = N2, y = result004, SE = 0, Algorithm = rep("P=20",length(N))  )
      result5 = data.frame( x = N2, y = result005, SE = 0, Algorithm = rep("P=25",length(N))  )
    }
    result = as.data.frame( rbind(result1, result2,result3, result4, result5) )
    
    result[,4] <- factor(result[,4], levels = c("P=5","P=10", "P=15", "P=20","P=25"))
    
    colnames(result)[4] = "Nodes"
    linesize= 1.0
    id = c(1,2,3,4,5,6,7)
    #id = c(1,2, 4:7)
    Color = c("black", "red", "blue", "purple4", "green4", "tomato4", "orange")[id]
    LineType = c("solid", "solid", "dashed", "dashed",  "longdash" , "longdash",  "dotted")[id]
    LineShape = c(16, 17, 18, 6, 7, 8, 12)[id]
    
    ### png 450/300 ###
    q = ggplot( data = result, aes(y= y, x= x, group = Nodes, color = Nodes, shape = Nodes) )
    #q = q + geom_errorbar(data = result, aes(ymin = y-SE ,  ymax = y+SE), width= 10)
    q = q + geom_line( aes(linetype=Nodes), size = linesize ) +  
      scale_color_manual(values=Color  ) +
      scale_linetype_manual(values= LineType) +
      scale_shape_manual(values = LineShape)
    q = q + geom_point( size = 2.5 ) 
    if(scale_param){
      q = q + xlab("C = h/log(p)") + ylab("") + ylim(0,1)
    }else{
      q = q + xlab("Trimmed Sample Size") + ylab("") + ylim(0,1)
    }
    #  
    q = q + theme_minimal(base_size = 16, base_family = "",
                          base_line_size = 11/22, base_rect_size = 10/100)
    q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
    
    print(q) 
    ############ Plot End #############
  
    dev.off()
  }
}

############# LTS Performce by H_ratio ###########
N = seq(50, 500, by = 50); d = 1
scale_param = FALSE
for(model in c("OutlierAll")){
  for(P in c(20)){
    for(num_out in c(1)){
    model_dir = file.path(result_root, model)
    if (!dir.exists(model_dir)) {
      warning(sprintf("Result directory '%s' not found; skipping.", model_dir))
      next
    }
    plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"P",P,"D",d,"B",num_out,"H.png")
    png(plot_filename, width = 450, height = 300)

    ################### plotting ############
    result1 = result2 = result3 = result4 = result5 = result6 = result7 = NULL
    result001 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result002 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.6, thresh = 2, result_dir = model_dir)[[5]]/100
    result003 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.7, thresh = 2, result_dir = model_dir)[[5]]/100
    result004 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.8, thresh = 2, result_dir = model_dir)[[5]]/100
    result005 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.9, thresh = 2, result_dir = model_dir)[[5]]/100
    
    if( scale_param ){
      result1 = data.frame( x = N*0.5, y = result001, SE = 0, Algorithm = rep("alpha=0.5",length(N))  )
      result2 = data.frame( x = N*0.6, y = result002, SE = 0, Algorithm = rep("alpha=0.6",length(N))  )
      result3 = data.frame( x = N*0.7, y = result003, SE = 0, Algorithm = rep("alpha=0.7",length(N))  )
      result4 = data.frame( x = N*0.8, y = result004, SE = 0, Algorithm = rep("alpha=0.8",length(N))  )
      result5 = data.frame( x = N*0.9, y = result005, SE = 0, Algorithm = rep("alpha=0.9",length(N))  )
    }else{
      result1 = data.frame( x = N, y = result001, SE = 0, Algorithm = rep("alpha=0.5",length(N))  )
      result2 = data.frame( x = N, y = result002, SE = 0, Algorithm = rep("alpha=0.6",length(N))  )
      result3 = data.frame( x = N, y = result003, SE = 0, Algorithm = rep("alpha=0.7",length(N))  )
      result4 = data.frame( x = N, y = result004, SE = 0, Algorithm = rep("alpha=0.8",length(N))  )
      result5 = data.frame( x = N, y = result005, SE = 0, Algorithm = rep("alpha=0.9",length(N))  )
    }
    result = as.data.frame(rbind(result1,result2,result3,result4,result5))
    
    result[,4] <- factor(result[,4], levels = c("alpha=0.5","alpha=0.6", "alpha=0.7", "alpha=0.8","alpha=0.9"))
    
    colnames(result)[4] = "Nodes"
    linesize= 1.0
    id = c(1,2,3,4,5,6,7)
    #id = c(1,2, 4:7)
    Color = c("black", "red", "blue", "purple4", "green4", "tomato4", "orange")[id]
    LineType = c("solid", "solid", "dashed", "dashed",  "longdash" , "longdash",  "dotted")[id]
    LineShape = c(16, 17, 18, 6, 7, 8, 12)[id]
    
    ### png 450/300 ###
    q = ggplot( data = result, aes(y= y, x= x, group = Nodes, color = Nodes, shape = Nodes) )
    #q = q + geom_errorbar(data = result, aes(ymin = y-SE ,  ymax = y+SE), width= 10)
    q = q + geom_line( aes(linetype=Nodes), size = linesize ) +  
      scale_color_manual(values=Color  ) +
      scale_linetype_manual(values= LineType) +
      scale_shape_manual(values = LineShape)
    q = q + geom_point( size = 2.5 ) 
    q = q + xlab("Sample Size") + ylab("") + ylim(0,1)
    
    #q = q + xlab("Sample Size") + ylab("") + ylim(0,1)
    q = q + theme_minimal(base_size = 16, base_family = "",
                          base_line_size = 11/22, base_rect_size = 10/100)
    q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
    
    print(q) 
    ############ Plot End #############
    
    dev.off()
  }
}}


############# LTS Performce by eta ###########
N = seq(50, 500, by = 50)
N2 = N/2
"Outlier1"

scale_param = FALSE
for(model in c("OutlierAll")){
  for(P in c(20)){
    for(num_out in c(1)){
      model_dir = file.path(result_root, model)
      if (!dir.exists(model_dir)) {
        warning(sprintf("Result directory '%s' not found; skipping.", model_dir))
        next
      }

      if(scale_param){
        plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"P",P,"D",d,"B",num_out,"EtaC.png")
      }else{
        plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"P",P,"D",d,"B",num_out,"Eta.png")
      }
      
      
      png(plot_filename, width = 450, height = 300)
      
      ################### plotting ############
      result1 = result2 = result3 = result4 = result5 = result6 = result7 = NULL
      result001 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 1, result_dir = model_dir)[[5]]/100
      result002 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 1.5, result_dir = model_dir)[[5]]/100
      result003 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
      result004 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 3, result_dir = model_dir)[[5]]/100
      result005 = simulation_data_fun(N, P, d = 1, "LTS", num_out = num_out, h_ratio = 0.5, thresh = 1000, result_dir = model_dir)[[5]]/100
    
      if( scale_param ){ # variance를 몰라서...
        result1 = data.frame( x = (N2- num_out) /pnorm(1) , y = result001, SE = 0, Algorithm = rep("eta=1.0",length(N))  )
        result2 = data.frame( x = (N2- num_out) /pnorm(1.5), y = result002, SE = 0, Algorithm = rep("eta=1.5",length(N))  )
        result3 = data.frame( x = (N2- num_out) /pnorm(2), y = result003, SE = 0, Algorithm = rep("eta=2.0",length(N))  )
        result4 = data.frame( x = (N2- num_out) /pnorm(3), y = result004, SE = 0, Algorithm = rep("eta=2.5",length(N))  )
        result5 = data.frame( x = (N2- num_out) /pnorm(100), y = result005, SE = 0, Algorithm = rep("eta=3.0",length(N))  )
      }else{
        result1 = data.frame( x = N2, y = result001, SE = 0, Algorithm = rep("eta=1",length(N))  )
        result2 = data.frame( x = N2, y = result002, SE = 0, Algorithm = rep("eta=1.5",length(N))  )
        result3 = data.frame( x = N2, y = result003, SE = 0, Algorithm = rep("eta=2",length(N))  )
        result4 = data.frame( x = N2, y = result004, SE = 0, Algorithm = rep("eta=3",length(N))  )
        result5 = data.frame( x = N2, y = result005, SE = 0, Algorithm = rep("eta=100",length(N))  )
      }
      
      result = as.data.frame(rbind(result1,result2,result3,result4,result5))
      
      result[,4] <- factor(result[,4], levels = c("eta=1","eta=1.5","eta=2", "eta=3","eta=100"))
      
      colnames(result)[4] = "Nodes"
      linesize= 1.0
      
      #id = c(1,2, 4:7)
      Color = c("black", "red", "blue", "purple4", "green4", "tomato4", "orange")
      LineType = c("dashed", "dashed", "solid", "longdash",  "longdash" , "longdash",  "dotted")
      LineShape = c(16, 17, 11, 15, 7, 8, 12)
      
      ### png 450/300 ###
      q = ggplot( data = result, aes(y= y, x= x, group = Nodes, color = Nodes, shape = Nodes) )
      #q = q + geom_errorbar(data = result, aes(ymin = y-SE ,  ymax = y+SE), width= 10)
      q = q + geom_line( aes(linetype=Nodes), size = linesize ) +  
        scale_color_manual(values=Color  ) +
        scale_linetype_manual(values= LineType) +
        scale_shape_manual(values = LineShape)
      q = q + geom_point( size = 2.5 ) 
      if(scale_param){
        q = q + xlab("C = n/log(p)") + ylab("") + ylim(0,1)
      }else{
        q = q + xlab("Trimmed Sample Size") + ylab("") + ylim(0,1)
      }
      #q = q + xlab("Sample Size") + ylab("") + ylim(0,1)
      q = q + theme_minimal(base_size = 16, base_family = "",
                            base_line_size = 11/22, base_rect_size = 10/100)
      q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
      
      print(q)
      ############ Plot End #############
      
      dev.off()
    }
  }}


############# LTS Performce by B ###########
N = seq(50, 500, by = 50); d = 1
scale_param = F
for(model in c("Outlier1")){
  for(P in c(20)){
    model_dir = file.path(result_root, model)
    if (!dir.exists(model_dir)) {
      warning(sprintf("Result directory '%s' not found; skipping.", model_dir))
      next
    }
    if(scale_param){
      plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"P",P,"D",d,"BC.png")
    }else{
      plot_filename = paste0("D:/OneDrive - UOS/2023Research/2022GaussianDAG_LTS/report/figures/",model,"P",P,"D",d,"B.png")
    }
      
    png(plot_filename, width = 450, height = 300)
    
    ################### plotting ############
    result1 = result2 = result3 = result4 = result5 = result6 = result7 = NULL
    result001 = simulation_data_fun(N, P, d = 1, "LTS", num_out = 30, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result002 = simulation_data_fun(N, P, d = 1, "LTS", num_out = 60, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result003 = simulation_data_fun(N, P, d = 1, "LTS", num_out = 90, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result004 = simulation_data_fun(N, P, d = 1, "LTS", num_out = 120, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    result005 = simulation_data_fun(N, P, d = 1, "LTS", num_out = 150, h_ratio = 0.5, thresh = 2, result_dir = model_dir)[[5]]/100
    
    if( scale_param ){
      result1 = data.frame( x = N/30, y = result001, SE = 0, Algorithm = rep("|B|=30",length(N))  )
      result2 = data.frame( x = N/60, y = result002, SE = 0, Algorithm = rep("|B|=60",length(N))  )
      result3 = data.frame( x = N/90, y = result003, SE = 0, Algorithm = rep("|B|=90",length(N))  )
      result4 = data.frame( x = N/120, y = result004, SE = 0, Algorithm = rep("|B|=120",length(N))  )
      result5 = data.frame( x = N/150, y = result005, SE = 0, Algorithm = rep("|B|=150",length(N))  )
    }else{
      result1 = data.frame( x = N, y = result001, SE = 0, Algorithm = rep("|B|=30",length(N))  )
      result2 = data.frame( x = N, y = result002, SE = 0, Algorithm = rep("|B|=60",length(N))  )
      result3 = data.frame( x = N, y = result003, SE = 0, Algorithm = rep("|B|=90",length(N))  )
      result4 = data.frame( x = N, y = result004, SE = 0, Algorithm = rep("|B|=120",length(N))  )
      result5 = data.frame( x = N, y = result005, SE = 0, Algorithm = rep("|B|=150",length(N))  )
    }
    result = as.data.frame(rbind(result1,result2,result3,result4,result5))
    result[,4] <- factor(result[,4], levels = c("|B|=30","|B|=60", "|B|=90", "|B|=120","|B|=150"))
    
    colnames(result)[4] = "Nodes"
    linesize= 1.0
    id = c(1,2,3,4,5,6,7)
    #id = c(1,2, 4:7)
    Color = c("black", "red", "blue", "purple4", "green4", "tomato4", "orange")[id]
    LineType = c("solid", "solid", "dashed", "dashed",  "longdash" , "longdash",  "dotted")[id]
    LineShape = c(16, 17, 18, 6, 7, 8, 12)[id]
    
    ### png 450/300 ###
    q = ggplot( data = result, aes(y= y, x= x, group = Nodes, color = Nodes, shape = Nodes) )
    #q = q + geom_errorbar(data = result, aes(ymin = y-SE ,  ymax = y+SE), width= 10)
    q = q + geom_line( aes(linetype=Nodes), size = linesize ) +  
      scale_color_manual(values=Color  ) +
      scale_linetype_manual(values= LineType) +
      scale_shape_manual(values = LineShape)
    q = q + geom_point( size = 2.5 ) 
    if(scale_param){
      q = q + xlab("C = h/|B|") + ylab("") + ylim(0,1)
    }else{
      q = q + xlab("Trimmed Sample Size") + ylab("") + ylim(0,1)
    }
    q = q + theme_minimal(base_size = 16, base_family = "",
                          base_line_size = 11/22, base_rect_size = 10/100)
    q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
    
    print(q) 
    ############ Plot End #############
    
    dev.off()
  }}



############## Run Time Result #################
setwd("D:/2020HighLSEM/Result/GSEM")

N = seq(100, 1000, by = 100)
P = c(100); d = 8


N = 1000
P = c(25, 50, 100, 150); d = 8
j = 6; r_index = 1


#result1[5,2] = 252.3896

###################
result1 = result2 = result3 = result4 =  result5 = result6 =NULL
result1 = simulation_data_fun(N, P, d, "HLSM", j)[[r_index]]
result4 = simulation_data_fun(N, P, d, "LISTEN", j)[[r_index]]
result5 = simulation_data_fun(N, P, d, "TD", j)[[r_index]]

result1[,1] = P
result4[,1] = P
result5[,1] = P
result = as.data.frame( rbind( result1, result4, result5 ) )
#colnames(result)[4] = "NodeSize"

linesize= 1.0
#id = 1:5
#id = c(1,2,3,4,5,6)
id = c(1,2,3,4,5,6,7)
Color = c("black", "red", "blue", "purple4", "green4", "tomato4", "orange")[id]
LineType = c("solid", "longdash", "dashed", "dotted",  "longdash" , "longdash",  "dotted")[id]
LineShape = c(16, 17, 18, 6, 7, 8, 12)[id]

### png 450/300 ###
q = ggplot( data = result, aes(y= log(y), x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
q = q + geom_errorbar(data = result, aes(ymin = log( y-SE ) ,  ymax = log( y+ SE) ), width= 10)
q = q + geom_line( aes(linetype=Algorithm), size = linesize ) +
  scale_color_manual(values=Color  ) +
  scale_linetype_manual(values= LineType) +
  scale_shape_manual(values = LineShape)
q = q + geom_point( size = 2.5 )
q = q + xlab("Node Size") + ylab("Log of Run Time") + ylim(1,8.5)
q = q + theme_minimal(base_size = 16, base_family = "",
                      base_line_size = 11/22, base_rect_size = 9)
q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "right" )
q

#############################




library(xtable)
result = cbind(result1[,1:3], result[,2:3])
xtable( (result2) )
x =  result1$y
y =  result1$y
z = result1$y
data = rbind( cbind(x, 20), cbind(y, 50), cbind(z, 80) )
mean(y)/mean(x)
/, mean(z) )
lm(data[,1] ~ data[,2])

linesize= 1.0
#id = 1:5
id = c(1,2,3,4)
Color = c("black", "red", "blue", "purple4", "green4", "tomato4")[id]
LineType = c("solid", "dotted", "dashed", "dashed",  "dashed" , "dashed")[id]
LineShape = c(16, 17, 18, 6, 7, 8)[id]

### png 450/300 ###
q = ggplot( data = result, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
q = q + geom_errorbar(data = result, aes(ymin = y-SE ,  ymax = y+SE), width= 10)
q = q + geom_line( aes(linetype=Algorithm), size = linesize ) +
  scale_color_manual(values=Color  ) +
  scale_linetype_manual(values= LineType) +
  scale_shape_manual(values = LineShape)
q = q + geom_point( size = 2.5 )
q = q + xlab("Sample Size") + ylab("")
q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
q
########################




#######################
MRS_result1 =MRS_result2 = Oracle_result = GES_result = MMHC_result = ODS_result1= NULL
MRS_result1 = simulation_data_fun(N, P, "MRS", 1)
#MRS_result2 = simulation_data_fun(N, P, "MRS", 1)
ODS_result1 = simulation_data_fun(N, P, "ODS", 1)
Oracle_result = simulation_data_fun(N, P, "MRS", 3)
GES_result = simulation_data_fun(N, P, "GES", 1)
MMHC_result = simulation_data_fun(N, P, "MMHC", 1)

#Oracle_result[[1]][,2] = 1
#Oracle_result[[1]][,3] = 0
#Oracle_result[[2]][,2] = 1
#Oracle_result[[2]][,3] = 0

MRS_result1[[1]][,4] = "2.MRS DAG"
MRS_result1[[2]][,4] = "2.MRS DAG"
ODS_result1[[1]][,4] = "3.ODS DAG"
ODS_result1[[2]][,4] = "3.ODS DAG"
Oracle_result[[1]][,4] ="1.Oracle"
Oracle_result[[2]][,4] ="1.Oracle"
GES_result[[1]][,4] ="4.GES DAG"
GES_result[[2]][,4] ="4.GES DAG"
MMHC_result[[1]][,4] ="5.MMHC DAG"
MMHC_result[[2]][,4] ="5.MMHC DAG"

MRS_result1[[1]][,4] = "2.MRS MEC"
MRS_result1[[2]][,4] = "2.MRS MEC"
ODS_result1[[1]][,4] = "3.ODS MEC"
ODS_result1[[2]][,4] = "3.ODS MEC"
Oracle_result[[1]][,4] ="1.Oracle"
Oracle_result[[2]][,4] ="1.Oracle"
GES_result[[1]][,4] ="4.GES MEC"
GES_result[[2]][,4] ="4.GES MEC"
MMHC_result[[1]][,4] ="5.MMHC MEC"
MMHC_result[[2]][,4] ="5.MMHC MEC"

#MRS_result2[[1]][,4] = "4.MRS MEC"
#MRS_result2[[2]][,4] = "4.MRS MEC"

result1 = rbind( Oracle_result[[1]], MRS_result1[[1]], MRS_result2[[1]], ODS_result1[[1]], GES_result[[1]], MMHC_result[[1]] )
result2 = rbind( Oracle_result[[1]], MRS_result1[[2]], MRS_result2[[2]], ODS_result1[[2]], GES_result[[2]], MMHC_result[[2]] )

######################
#result1 = rbind( Oracle_result[[1]], MRS_result1[[1]], MRS_result2[[1]], GES_result[[1]], MMHC_result[[1]] )
#result2 = rbind( Oracle_result[[1]], MRS_result1[[2]], MRS_result2[[2]], GES_result[[2]], MMHC_result[[2]] )

linesize= 1.0
id = 1:5
Color = c("black", "red", "blue", "purple4", "green4", "tomato4")[id]
LineType = c("solid", "solid", "solid", "dashed",  "dashed" , "dashed")[id]
LineShape = c(16, 17, 18, 6, 7, 8)[id]
### png 450/300 ###
q = ggplot( data = result1, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
#q = ggplot( data = result_time, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
q = q + ylim(0.0, 1.00)
q = q + xlab("Sample Size") + ylab("")
q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "right" )
q

q = ggplot( data = result2, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
#q = ggplot( data = result_time, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
q = q + ylim(0.0, 1.00)
q = q + geom_errorbar(data = result2, aes(ymin = y-SE ,  ymax = y+SE), width= 10)
q = q + geom_line( aes(linetype=Algorithm), size = linesize ) +
  scale_color_manual(values=Color  ) +
  scale_linetype_manual(values=LineType) +
  scale_shape_manual(values = LineShape)
q = q + geom_point( size = 2.5 )
q = q + xlab("Sample Size") + ylab("")
q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
q

### png 450/300 ###

q = ggplot( data = result1, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
#q = ggplot( data = result_time, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
q = q + ylim(0.0, 1)
q = q + geom_errorbar(data = result1, aes(ymin = y-SE ,  ymax = y+SE), width= 40)
q = q + geom_line( aes(linetype=Algorithm), size = linesize ) +
  scale_color_manual(values=c("black", "red", "purple4", "green4", "tomato4")  ) +
  scale_linetype_manual(values=c("solid", "solid", "dashed",  "dashed" , "dashed")) +
  scale_shape_manual(values = c(16, 17, 6, 7, 8))
q = q + geom_point( size = 2.5 )
q = q + xlab("Sample Size") + ylab("")
q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
q

### png 450/300 ###

q = ggplot( data = result2, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
#q = ggplot( data = result_time, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
q = q + ylim(0.0, 1)
q = q + geom_errorbar(data = result2, aes(ymin = y-SE ,  ymax = y+SE), width= 40)
q = q + geom_line( aes(linetype=Algorithm), size = linesize ) +
  scale_color_manual(values=c("black", "red", "purple4", "green4", "tomato4")  ) +
  scale_linetype_manual(values=c("solid", "solid", "dashed",  "dashed" , "dashed")) +
  scale_shape_manual(values = c(16, 17, 6, 7, 8))
q = q + geom_point( size = 2.5 )
q = q + xlab("Sample Size") + ylab("")
q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
q

q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "right" )
q


###


setwd("C:/2018AISTATS_GSEM/R-code/PoissonDAGresult_logLink_d1/")
setwd("C:/2018AISTATS_GSEM/R-code/PoissonDAGresult_logLink_d2/")
setwd("C:/2018AISTATS_GSEM/R-code/PoissonDAGresult_logLink_d5/")
setwd("C:/2018AISTATS_GSEM/R-code/PoissonDAGresult_logLink_d10/")


N = 500
P = seq(10, 200, by = 10)

N = seq(100,1000, by=100)
P = c(500)

MRS_time = simulation_data_fun(N, P, "MRS", 6)
GES_time = simulation_data_fun(N, P, "GES", 6)
MMHC_time = simulation_data_fun(N, P, "MMHC", 6)

MRS_time[[1]][,4] = "1.MRS"
GES_time[[1]][,4] ="2.GES"
MMHC_time[[1]][,4] ="3.MMHC"

MRS_time[[1]]$x = P
GES_time[[1]]$x = P
MMHC_time[[1]]$x = P

Oracle_time$y = 2*log(P) - 2*log(10)
Oracle_time[,4] = "4.2LogP"
Oracle_time2 = MRS_time[[1]]
Oracle_time2$y = 3*log(P) - 3*log(10)
Oracle_time2[,4] = "5.2LogP"

result = rbind( MRS_time[[1]], GES_time[[1]], MMHC_time[[1]])
result$y = (result$y)
#result$SE = 0
q = ggplot( data = result, aes(y= y, x= x, group = Algorithm, color = Algorithm, shape = Algorithm) )
q = q + geom_errorbar(data = result, aes(ymin = y-SE ,  ymax = y+SE), width= 50)
q = q + geom_line( aes(linetype=Algorithm), size = 1 ) +
  scale_color_manual(values=c("red", "green4", "tomato4", "blue", "purple")  ) +
  scale_linetype_manual(values=c("solid", "dashed",  "dashed",  "dashed",  "dashed" )) +
  scale_shape_manual(values = c(17, 7, 8, 9, 10))
q = q + geom_point(size = 2.5)
q = q + xlab("Sample Size") + ylab("Run Time (Seconds)") #+ ggtitle("Skeleton")
q = q + theme(title= element_text(face = "bold", color= "black"), legend.position = "none" )
q



