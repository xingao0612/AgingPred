#' AgingPred
#'
#' @param text_size Text size
#'
#' @return A shiny app
#'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
#' @examples
#' \dontrun{
#' AgingPred()
#' }
## Author: Xin Gao
## Maintainer: Xin Gao <gaoxin_0612@163.com>
## Description: A Shiny App For Aging Prediction.
### All rights reserved. Please contact the author for permission before using this code for any purpose.

AgingPred <- function(text_size = 20) {
  options(shiny.maxRequestSize = 550*1024^2) ###限制上传表达谱内存大小
  #options(repos = c(BiocManager::repositories(), CRAN = "https://cloud.r-project.org/")) ####不然无法部署limma这类包

  # 创建UI界面
  ui <- dashboardPage(
    dashboardHeader(title = "Aging Prediction"),
    dashboardSidebar(
      fileInput("expression_file", "Upload Expression Profile File (in CSV format)"),

      fileInput("clin_file", "Upload Sample Information File (in CSV format)"),

      actionButton("calculate_button", "Start Prediction")
    ),

    dashboardBody(
      fluidRow(
        box(
          title = "Correlation between Predicted Age and Actual Age", status = "primary", solidHeader = TRUE,
          collapsible = TRUE,  withSpinner(plotOutput("age_plot", height = 350), image = "loading.gif"),
          downloadButton("age_plot1", "Download Plot"),
          downloadButton("downloadData1", "Download Data")),


        box(
          title = "Stacked Bar Plot of Aging Type Distribution", status = "primary", solidHeader = TRUE,
          collapsible = TRUE, withSpinner(plotOutput("bar_plot", height = 350), image = "loading.gif"),
          downloadButton("bar_plot1", "Download Plot"),
          downloadButton("downloadData2", "Download Data"))

      ),


      fluidRow(
        box(status = "success",
            h3(tagList(shiny::icon("question-circle"), strong("About this software"))),
            h5(
              "Welcome to", strong("AgingPred!"), "AgingPred includes random forest prediction models for age and aging rate, which can be used for age and aging rate analysis in both healthy individuals and disease populations, particularly patients with infectious diseases. (Due to potential biases in model training, the predicted results are for scientific research reference only and are not yet suitable for clinical applications.)", br(), br(),
              strong("Affiliation: "), "Peking Union Medical College, Chinese Academy of Medical Sciences; The Key Laboratory of Geriatrics, Beijing Institute of Geriatrics, National Health Commission, Beijing Hospital. (中国医学科学院&北京协和医学院；北京医院，国家卫生健康委员会北京老年医学研究所老年医学重点实验室)", br(),
              strong("Author: "), "Xin Gao (高鑫), E-mail: ", span("gaoxin_0612@163.com", style = "color:blue"), br(),
              "Please cite XXXXX articles after using this software.", br(), br(),
              strong("Upload data description:"),
              "The uploaded expression profile must include the genes ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=CD248", "CD248"),
              ", ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=PHYKPL", "PHGDH"),
              ", ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=PHYKPL", "PHYKPL"),
              ", ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=BPGM", "BPGM"),
              ", ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=SFXN2", "SFXN2"),
              ", ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=TWF1", "TWF1"),
              ", ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=MXRA8", "MXRA8"),
              ", ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=NOG", "NOG"),
              ", ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=TTC24", "TTC24"),
              ", and ",
              tags$a(href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=CACHD1", "CACHD1"),
              ", which are part of the random forest prediction model. Please click on the gene name to view the gene aliases!",
              br(),
              "It is recommended to use RNA-seq data from human peripheral blood leukocytes. The required example data can be downloaded from the link below. If the sample information file does not contain any grouping information, please fill the 'status' column with 'none'.",
              style = "text-align:left;color:black;background-color:lavender;padding:15px;border-radius:10px"
            ),
            actionButton("openButton" , label = HTML(paste0(icon("cloud-download"), "Example Expression Profile File"))),
            actionButton("openButton2" , label = HTML(paste0(icon("cloud-download"), "Example Sample Information File"))),
            width = 12
        )
      )



    )







  )



  # 创建Server逻辑#################################################################################################
  server <- function(input, output, session) {

    ###提供示范数据下载
    observeEvent(input$openButton, {
      browseURL("https://figshare.com/ndownloader/files/40473161") #表达谱
    })

    observeEvent(input$openButton2, {
      browseURL("https://figshare.com/ndownloader/files/40473170") #年龄分组信息
    })

    # 目标基因列表
    target_genes <- c("CD248", "PHGDH", "PHYKPL", "BPGM", "SFXN2", "TWF1", "MXRA8", "NOG"
                      , "TTC24", "CACHD1")
    # 处理上传的表达谱文件
    rt <- reactive({
      req(input$expression_file)
      read.csv(input$expression_file$datapath, row.names = 1)

    })

    # 检查是否包含目标基因
    check_target_genes <- function(rt) {
      has_target_genes <- all(target_genes %in% rownames(rt))
      return(has_target_genes)
    }

    # 读取临床信息
    clin <- reactive({
      req(input$clin_file)
      read.csv(input$clin_file$datapath)

    })

    # 响应点击按钮事件，计算年龄并绘制散点图
    observeEvent(input$calculate_button, {
      if (!check_target_genes(rt())) {
        showModal(modalDialog(
          title = "Error",
          "The target genes (CD248, PHGDH, HYKPL, BPGM, SFXN2, TWF1, MXRA8, NOG, TTC24, CACHD1) are missing in the expression profile! Please check the uploaded expression profile file.",
          easyClose = TRUE
        ))
      } else {
        ###########主要执行程序###############################################################
        library(limma)
        library(sva)
        #rt=read.csv("geneMatrix.csv",row.names = 1) ##上传数据
        usersamples=colnames(rt())
        load(system.file("illumina.merge.combat.Rdata", package = "AgingPred"))
        intersectGenes=intersect(rownames(rt()), rownames(rt_my))

        ####对上传的数据预处理###################
        rt=as.matrix(rt())

        #####开发模式############################################################
        #rt=read.csv("geneMatrix.csv", row.names = 1)
        #rt=as.matrix(rt)
        #usersamples=colnames(rt)
        #load("illumina.merge.combat.Rdata")
        #intersectGenes=intersect(rownames(rt), rownames(rt_my))
        #clin1=read.csv("clinical.csv")
        ###############################################################################

        exp=rt[,1:ncol(rt)]
        dimnames=list(rownames(exp),colnames(exp))
        data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
        rt=avereps(data)
        colnames(rt)=paste0(colnames(rt))

        #对数值大的数据取log2
        qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
        if(LogC){
          rt[rt<0]=0
          rt=log2(rt+1)}

        rt=normalizeBetweenArrays(rt)
        rt=as.data.frame(rt)
        #数据合并
        rt=rt[intersectGenes,]
        rt_my=rt_my[intersectGenes,]

        allTab=cbind(rt,rt_my)
        ###构建批次信息
        batchType=c()
        batchType=c(batchType, rep(1,ncol(rt)))
        batchType=c(batchType, rep(2,ncol(rt_my)))

        #对数据进行矫正，输出矫正后的结果
        outTab=ComBat(allTab, batchType, par.prior=TRUE)

        ##提取单个数据集矫正后的数据

        rt=as.data.frame(outTab[,usersamples])

        ##########验证模型3####################################################################
        #randomForest 包的随机森林
        library(randomForest)
        library(ggplot2)
        library(ggpubr)
        ###加载内置数据
        load(system.file("otu_train.select.forest.Rdata", package = "AgingPred"))
        load(system.file("training.pred.result.Rdata", package = "AgingPred"))
        load( system.file("test1.pred.result.Rdata", package = "AgingPred"))

        #rt=read.csv("GSE185576.combat.csv", row.names = 1)
        rt=as.data.frame(t(rt))
        rt$sample=rownames(rt)

        clin1=clin()  ###上传临床信息
        clin2=clin1[,c(1,2)]
        clin2=as.data.frame(clin2)
        clin2$age=as.numeric(clin2$age)

        dat3=merge(rt,clin2,by.x = "sample", by.y = "sample")
        rownames(dat3)=dat3[,1]

        dat3=dat3[,-1]

        testdata3=dat3[,1:c(ncol(dat3)-1)]

        #使用测试集，查看预测精度
        predict_test <- predict(otu_train.select.forest, testdata3)

        c=cbind(act_age=dat3$age,predict.age=predict_test, delta_age=c(predict_test-dat3$age))
        c=as.data.frame(c)

        df3=as.data.frame(cbind(dat3$age, predict_test))
        df3$sample=rownames(df3)
        df3=merge(df3, clin1, by.x = "sample", by.y = "sample")
        head(df3)

        p1 <-ggplot(df3, aes(age, predict_test,color = status)) +
          xlab("Actual age (years)")+ylab("Predicted age (years)")+
          geom_point(aes(color = status), alpha=0.8)+ geom_smooth(method="lm",formula = y ~ x, color="#333333", fill = "#cbc9e2") + theme_test()+
          stat_cor(method = 'spearman', aes(x =age, y =predict_test), show.legend = FALSE)+geom_text(data = df3,
                                                                                                     aes(label = sprintf("R = %.2f\nP = %.3f",
                                                                                                                         cor(age, predict_test, method = "spearman"),
                                                                                                                         cor.test(age, predict_test, method = "spearman")$p.value)),
                                                                                                     x = max(df3$age), y = max(df3$predict_test), hjust = 1, vjust = 1, color = "#666666")
        +
          scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"))


        output$age_plot <- renderPlot({
          p1
        })
        ####根据最相关文献的方法计算误差矫正delta年龄###################################################################3
        library(dplyr)
        a=as.data.frame(a)
        b=as.data.frame(b)
        c=as.data.frame(c)
        a$group=0
        b$group=1
        c$group=2

        abc=rbind(a,b,c)

        e=loess(delta_age~act_age+group, abc)
        f=predict(e, newdata = abc)
        curated_delta_age=abc$delta_age-f
        curated_delta_age=as.data.frame(curated_delta_age)

        all_pred=cbind(abc, curated_delta_age)

        #################################################################################################
        #####NGDC样本单独划分四分位数
        ngdc=subset(all_pred, group=="2")
        del=which(colnames(ngdc)=="group")
        ngdc=ngdc[,-del]###去掉运行标记列

        quantile(ngdc$curated_delta_age)
        q1=quantile(ngdc$curated_delta_age, probs = 0.25)
        q2=quantile(ngdc$curated_delta_age, probs = 0.50)
        q3=quantile(ngdc$curated_delta_age, probs = 0.75)

        ngdc$Aging_group <- case_when(ngdc$curated_delta_age > q3 ~ "quick-aging",
                                      ngdc$curated_delta_age < q1 ~ "slow-aging",
                                      ngdc$curated_delta_age <= q3 | ngdc$curated_delta_age >=q3 ~ "average-aging")


        #####绘制堆积图
        ngdc$sample=rownames(ngdc)
        dat=merge(ngdc, clin1, by.x = "sample", by.y = "sample")

        ##先把数据变成2x2列联表，然后用 chisq.test函数做.
        mytable <- table(dat$Aging_group,dat$status)
        mytable
        chisq.test(mytable,correct = F)

        mytable2=as.data.frame(mytable)

        colnames(mytable2)[1]="Type"
        colnames(mytable2)[2]="status"
        # 计算每个样本中的基因表达量百分比
        plot_data <- mytable2 %>%
          group_by(status) %>%
          mutate(Percent = Freq / sum(Freq))


        # 绘制百分比堆积图
        p2 <-ggplot(plot_data, aes(x = status, y = Percent, fill = Type )) +
          geom_col(position = "fill")+
          scale_fill_manual(values = c("#8DA0CB", "#FC8D62", "#66C2A5"))+xlab("")+
          geom_text(aes(label = scales::percent(Percent)), position = position_fill(vjust = 0.5)) +  # 添加百分比数值
          theme_test()

        output$bar_plot <- renderPlot({
          p2
        })
        ##########数据和图像下载单元#######################################################################################################


        # 添加PDF下载按钮
        #下载图像
        output$age_plot1 <- downloadHandler(
          filename = function() {
            "age_cor_plot.pdf"
          },
          content = function(file) {
            pdf(file, width = 6, height = 4.5)
            print(p1)
            dev.off()
          }
        )

        output$bar_plot1 <- downloadHandler(
          filename = function() {
            "stacked_bar_chart.pdf"
          },
          content = function(file) {
            pdf(file, width = 6, height = 4)
            print(p2)
            dev.off()
          }
        )
        # 下载结果数据
        output$downloadData1 <- downloadHandler(
          filename = function() {
            "predict_results.csv"
          },
          content = function(file) {
            write.csv(dat, file, row.names = FALSE)
          }
        )

        output$downloadData2 <- downloadHandler(
          filename = function() {
            "aging_percent.csv"
          },
          content = function(file) {
            write.csv(plot_data, file, row.names = FALSE)
          }
        )


      }
    })
  }

  shinyApp(ui, server)
}
