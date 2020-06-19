# plot 1: enviromental variables panel
                    require(ggplot2)
                    ggplot() + 
                      scale_x_discrete(expand = c(0,0))+
                      scale_y_discrete(expand = c(0,0))+
                      scale_fill_gradientn(colours= rainbow(15))+
                      geom_tile(data=out,aes(x=reorder(env,Freq), y=reorder(env.variable,Freq),fill=Freq))+
                      ylab('Envirotype ID\n')+ 
                      xlab("\nEnvironment ID")+
                      labs(fill='Frequency')+
                      theme(axis.title = element_text(size=19),
                            legend.text = element_text(size=9),
                            strip.text.y  = element_text(size=13,angle=360),
                            legend.title = element_text(size=17),
                            strip.background = element_rect(fill="gray95",size=1),
                            legend.position = 'bottom')
                    
                    # plot 2: distribution of envirotypes
                    ggplot() + 
                      #theme_void()+
                      scale_x_discrete(expand = c(0,0))+
                      scale_y_continuous(expand = c(0,0))+
                      # theme_pubclean()+
                      #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
                      facet_grid(var~interval)+
                      geom_bar(data=out, aes(y=Freq, x=env,fill=env.variable), 
                               position = "fill",stat = "identity",width = 1,size=.2)+
                      # scale_y_continuous(labels = scales::percent,expand = c(0,0))+ #coord_flip()+
                      ylab('Absolute Frequency\n of Occurence\n')+ 
                      xlab("\nEnvironment ID")+
                      labs(fill='Envirotype')+
                      theme(axis.title = element_text(size=19),
                            #   axis.text = element_blank(),
                            legend.text = element_text(size=9),
                            strip.text = element_text(size=17),
                            legend.title = element_text(size=17),
                            strip.background = element_rect(fill="gray95",size=1),
                            legend.position = 'bottom')
