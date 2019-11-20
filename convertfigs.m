function convertfigs()
fprintf("Converting figures ...");
close all;
files = dir('figures');
for i = 1:length(files)
    currentfile = files(i).name;
    if(strcmp(currentfile, '.'))
        continue
    end
    if(strcmp(currentfile, '..'))
        continue
    end
    f = openfig(strcat('figures/',currentfile));
    %print(f,'-dpng',[currentfile(1:end-3),'png']);
    
    
    set(f,'PaperPosition',[0 0 2.4 1.4])
    set(findall(f,'-property','FontSize'),'FontSize',8)
    set(findall(f,'-property','LineWidth'),'LineWidth',1)
    outname = sprintf("eps/%seps", currentfile(1:end-3));
    print(gcf, outname, '-depsc','-r600')
    close(f);
end
fprintf("\n");
end