# plot 1: enviromental variables panel
require(ggplot2)
require(EnvRtype)

names.window <- c('1-intial growing','2-leaf expansion I','3-leaf expansion II','4-flowering','5-grain filling','6-maturation')
time.window  <- c(0,15,35,65,90,120)
out=EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T,time.window = time.window, names.window = names.window)

# plot 1
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

out=EnvTyping(env.data = df.clim,env.id = 'env',var.id='VPD',by.interval = T,time.window = time.window, names.window = names.window)
require(reshape2)
out=data.frame(out,colsplit(out$env.variable,pattern = '_',names = c('var2','freq','stages')))

ggplot() + 
  #theme_void()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # theme_pubclean()+
  #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
  facet_grid(var~interval)+ coord_flip()+
  geom_bar(data=out, aes(y=Freq, x=env,fill=freq), 
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


out=EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T,time.window = time.window, names.window = names.window)
require(reshape2)
out=data.frame(out,colsplit(out$env.variable,pattern = '_',names = c('var2','freq','stages')))

ggplot() + 
  #theme_void()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # theme_pubclean()+
  #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
  facet_grid(var~interval)+ coord_flip()+
  geom_bar(data=out, aes(y=Freq, x=env,fill=freq), 
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


out=EnvTyping(env.data = df.clim,env.id = 'env',var.id='SRAD',by.interval = T,time.window = time.window, names.window = names.window)
require(reshape2)
out=data.frame(out,colsplit(out$env.variable,pattern = '_',names = c('var2','freq','stages')))

ggplot() + 
  #theme_void()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # theme_pubclean()+
  #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
  facet_grid(var~interval)+ coord_flip()+
  geom_bar(data=out, aes(y=Freq, x=env,fill=freq), 
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


## another suggestion


ggplot() + 
  #theme_void()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # theme_pubclean()+
  #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
  #facet_grid(var~interval)+ 
  coord_flip()+
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


# plot 3: distribution of envirotypes per environment

out=EnvTyping(env.data = df.clim,env.id = 'env',var.id='VPD',by.interval = T,time.window = time.window, names.window = names.window)
require(reshape2)
out=data.frame(out,colsplit(out$env.variable,pattern = '_',names = c('var2','freq','stages')))

ggplot() + 
  #theme_void()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # theme_pubclean()+
  #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
  facet_grid(~env)+ #coord_flip()+
  geom_bar(data=out, aes(y=Freq, x=interval,fill=freq), 
           position = "fill",stat = "identity",width = 1,size=.2)+
  # scale_y_continuous(labels = scales::percent,expand = c(0,0))+ #coord_flip()+
  ylab('Absolute Frequency\n of Occurence\n')+ 
  xlab("\nEnvironment ID")+
  labs(fill='Envirotype')+
  theme(axis.title = element_text(size=19),
           axis.text.x  = element_text(hjust=1,angle=45),
        legend.text = element_text(size=9),
        strip.text = element_text(size=17),
        legend.title = element_text(size=17),
        strip.background = element_rect(fill="gray95",size=1),
        legend.position = 'bottom')


out=EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T,time.window = time.window, names.window = names.window)
require(reshape2)
out=data.frame(out,colsplit(out$env.variable,pattern = '_',names = c('var2','freq','stages')))

ggplot() + 
  #theme_void()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # theme_pubclean()+
  #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
  facet_grid(~env)+ #coord_flip()+
  geom_bar(data=out, aes(y=Freq, x=interval,fill=freq), 
           position = "fill",stat = "identity",width = 1,size=.2)+
  # scale_y_continuous(labels = scales::percent,expand = c(0,0))+ #coord_flip()+
  ylab('Absolute Frequency\n of Occurence\n')+ 
  xlab("\nEnvironment ID")+
  labs(fill='Envirotype')+
  theme(axis.title = element_text(size=19),
        axis.text.x  = element_text(hjust=1,angle=45),
        legend.text = element_text(size=9),
        strip.text = element_text(size=17),
        legend.title = element_text(size=17),
        strip.background = element_rect(fill="gray95",size=1),
        legend.position = 'bottom')

out=EnvTyping(env.data = df.clim,env.id = 'env',var.id='SRAD',by.interval = T,time.window = time.window, names.window = names.window)
require(reshape2)
out=data.frame(out,colsplit(out$env.variable,pattern = '_',names = c('var2','freq','stages')))

ggplot() + 
  #theme_void()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # theme_pubclean()+
  #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
  facet_grid(~env)+ #coord_flip()+
  geom_bar(data=out, aes(y=Freq, x=interval,fill=freq), 
           position = "fill",stat = "identity",width = 1,size=.2)+
  # scale_y_continuous(labels = scales::percent,expand = c(0,0))+ #coord_flip()+
  ylab('Absolute Frequency\n of Occurence\n')+ 
  xlab("\nEnvironment ID")+
  labs(fill='Envirotype')+
  theme(axis.title = element_text(size=19),
        axis.text.x  = element_text(hjust=1,angle=45),
        legend.text = element_text(size=9),
        strip.text = element_text(size=17),
        legend.title = element_text(size=17),
        strip.background = element_rect(fill="gray95",size=1),
        legend.position = 'bottom')
