import re
#import matplotlib as plt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as lines
import numpy as np
import imageio
import shutil
#data_path = "ANXA1.gff"
data_path = "GCF_000001405.39_GRCh38.p13_genomic.gtf"
f = open(data_path, 'r')
genetxt_path = './tmp.txt'
fp = open(genetxt_path, 'w')
allLines = f.readlines()
#genename = 'ANXA10'
#genename = 'FAM138A'
genename = 'ANKRD30A'
transcript_list = [] 
transcript_num = 0 
line_width = 5
for eachLine in allLines:
    if not eachLine.startswith('#'):
        eachLine_list = eachLine.split('\t')
        #name = re.search('gene_name "(.*?)";', eachLine_list[-1]).group(1)
        #name = re.search('gene=(.*?);', eachLine_list[-1])
        name = re.search('gene_id "(.*?)";', eachLine_list[-1])
        if name is not None:
            name = name.group(1)
        else:
            continue
        if name == genename:
            gene_exist = 1
            fp.write(eachLine)
            #if 'transcript_id' in eachLine:
                #import ipdb
                #ipdb.set_trace()
            if eachLine_list[2] == 'transcript':
                transcript_num += 1
                #transcript_name = re.search('transcript_name "(.*?)";', eachLine_list[-1]).group(1)
                transcript_name = re.search('transcript_id "(.*?)";', eachLine_list[-1])
                if transcript_name is not None:
                    transcript_name = transcript_name.group(1)
                else:
                    continue
                try:
                    #transcript_name = re.search('transcript_id=(.*?)\n', eachLine_list[-1]).group(1)
                    if transcript_name not in transcript_list:
                        transcript_list.append(transcript_name)
                except:
                    import ipdb
                    ipdb.set_trace()
f.close()
fp.close()
#若gtf文件中不存在相应基因的信息，弹出提示
if gene_exist == 0:
    self.failedFm = Frame(self.top)
    self.faillb = Label(self.failedFm, text="genename no exist,please check and retry").pack()
    self.failedFm.pack()
    self.top.update()
    time.sleep(3)
    self.genename.set('')
    self.failedFm.destroy()
    shutil.rmtree(result_path)
#反之向下执行
else:
    #import ipdb
    #ipdb.set_trace()
    fig = plt.figure(1)
    num = 0
    tmp_colors = ['lime', 'red', 'blue', 'yellow', 'yellow', 'w']
    names_tmp_colors = ['gene', 'CDS', 'exon', 'three_prime_utr', 'five_prime_utr', 'stop_codon']
    colors_legend_name = ['gene', 'CDS_exon', 'non_CDS_exon', 'UTR_exon']
    color_dict = dict(zip(names_tmp_colors, tmp_colors))
    result_path = './'
    png_path = result_path + '/' + genename + '.png'
    #读取之前保存的包含基因信息的txt文件，若不存在，弹出警告信息
    fp = open(genetxt_path, 'r')
    allLines = fp.readlines()
    for eachLine in allLines:
        #import ipdb
        #ipdb.set_trace()
        eachLine_list = eachLine.split('\t')
        if eachLine_list[2] == 'gene':
            #判断方向
            if eachLine_list[6] == '+':
                arr = '->'
            else:
                arr = '<-'

            ax = fig.add_axes([0.2, 0.2, 0.5, 0.6])
            arrow = mpatches.FancyArrowPatch(
                (int(eachLine_list[3]), 0.1),
                (int(eachLine_list[4]), 0.1),
                arrowstyle=arr,
                mutation_scale=25, lw=1, color='lime', antialiased=True)
            # 画箭头
            ax.add_patch(arrow)
            # 坐标轴标签
            ax.set_xlim(int(eachLine_list[3]), int(eachLine_list[4]))
            ax.set_ylim(-0.5, transcript_num + 1)
            ax.set_xticks(np.linspace(int(eachLine_list[3]), int(eachLine_list[4]), 5))
            ax.set_yticks([0.1] + list(range(1, transcript_num + 1)))
            #import ipdb
            #ipdb.set_trace()
            ax.set_yticklabels(['gene'] + transcript_list)
            #ax.set_yticklabels([0.1] + list(range(1, transcript_num + 1)))
            ax.set_xticklabels([str(j) for j in [int(i) for i in np.linspace (int(eachLine_list[3]), int(eachLine_list[4]),5)]])
            # 坐标轴显示
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            ax.get_xaxis().set_tick_params(direction='out')
            ax.tick_params(axis=u'y', which=u'both', length=0)  # 纵坐标刻度线不显示（length=0）
            # 坐标轴字体大小
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(6)


        elif eachLine_list[2] == 'transcript':
            num = num + 1
            line1 = [(int(eachLine_list[3]), num), (int(eachLine_list[4]), num)]
            (line1_xs, line1_ys) = zip(*line1)
            ax.add_line(lines.Line2D(line1_xs, line1_ys, linewidth=0.2,
                                    solid_capstyle='butt', solid_joinstyle='miter',
                                    antialiased=False, color='black'))

        elif eachLine_list[2] in color_dict.keys():
            # 添加结构图
            line2 = [(int(eachLine_list[3]) - 0.5, num), (int(eachLine_list[4]) + 0.5, num)]
            (line2_xs, line2_ys) = zip(*line2)
            ax.add_line(lines.Line2D(line2_xs, line2_ys,
                                    solid_capstyle='butt', solid_joinstyle='miter',
                                    linewidth=int(line_width), alpha=1,
                                    color=color_dict[eachLine_list[2]],
                                    antialiased=False))
        # 做图例
        # add_axes 是在一张图上指定特定区域作图，第一个数字为从左边%74处，下面20%处开始，宽20%，高60%区域作图
        ax_legend = fig.add_axes([0.76, 0.2, 0.2, 0.6])
        for i in range(len(colors_legend_name)):
            line3 = [(0, (9 - i) * 0.1), (0.1, (9 - i) * 0.1)]
            (line3_xs, line3_ys) = zip(*line3)
            ax_legend.add_line(lines.Line2D(line3_xs, line3_ys, linewidth=5,
                                            color=color_dict[names_tmp_colors[i]],
                                            solid_capstyle='butt', solid_joinstyle='miter',
                                            antialiased=False))
            ax_legend.text(0.2, (8.9 - i) * 0.1, colors_legend_name[i], fontsize=6)
            ax_legend.set_axis_off()
        #fig.suptitle('\n\n\nchr' + str(eachLine_list[0]) + ': ' + genename, fontsize=10)
        fig.suptitle('\n\n\n' + genename , fontsize=10)
        # 保存图片
        fig.savefig(png_path, dpi=150)
