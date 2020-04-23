


lib <- matrix(nrow = 1, ncol = 11)
colnames(lib) <- c("pos1","pos2","pos3","pos4","pos5","pos6","pos7","pos8","pos9","pos10","pos11")
lib[1, ] <- c(0,0,0,0,0,0,0,0,0,0,0)
lib <- as.data.frame(lib)

for (pos in 1:length(lib)) {
  for (possition in Mupexi_data_only_resonds_u_dup$peptide_position) {
    if (grepl(":", possition)==T) {
      pos1 <-  as.numeric(str_split(possition, ':')[[1]][1]  )
      pos2 <-  as.numeric( str_split(possition, ':')[[1]][2]  )
      posistions <- pos1:pos2
      if (pos %in% posistions) {
        lib[,pos] = lib[,pos]+1}
    }
    else if (grepl(":", possition)==F) {
      print(possition)
      possition <- as.numeric(possition)
      if (pos %in% possition) {
        lib[,pos] = lib[,pos]+1 }
    }
  }
  
}

lib <- as.data.frame(t(lib))
lib$pos <- rownames(lib)

ggplot(lib, aes(x=factor(pos,
                         levels = c("pos1" , "pos2" , "pos3" , "pos4" , "pos5" , "pos6" ,"pos7" , "pos8" , "pos9" , "pos10", "pos11")),
                y=V1)) +
  geom_bar(stat = "identity") + #stat = "identity", position = 'dodge', alpha = 0.5
  #  scale_fill_manual(values = fill_man)  +
  #  scale_y_log10(breaks = c(10, 50, 100,500, 1000, 5000, 10000), minor_breaks = NULL,
  #               labels = c('10', '50','100' ,'500','1000', '5000', '10000')) +
  labs(x = "", y ="number changes") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(size=18))