---
params: 
  set_title:
    label: "report title"
    value: "Project"
    input: text
  set_author: 
    label: "author"
    value: "Miles Smith"
    input: text
title: "Shakedown run"
author: "Miles Smith"
date: "29 September, 2021"
output: 
  html_document: 
    highlight: pygments
    toc: yes
    toc_float: yes
    keep_md: yes
  pdf_document:
    highlight: "pygments"
    toc: TRUE
    toc_depth: 3
always_allow_html: true
---

<style type="text/css">
    body .main-container {
      max-width: 1500px !important;
      width: 1500px !important;
    }
    body {
      max-width: 1500px !important;
    }
    caption {
      color: black;
      font-weight: bold;
      font-size: 1.0em;
    }
</style>





<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Disease_Class </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:left;"> percent </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;width: 10em; font-weight: bold;border-right:1px solid;"> control </td>
   <td style="text-align:right;width: 10em; border-right:1px solid;"> 87 </td>
   <td style="text-align:left;width: 10em; "> 23.77% </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 10em; font-weight: bold;border-right:1px solid;"> sle </td>
   <td style="text-align:right;width: 10em; border-right:1px solid;"> 279 </td>
   <td style="text-align:left;width: 10em; "> 76.23% </td>
  </tr>
</tbody>
</table>

\newpage

# Differential Gene Expression {.tabset .tabset-fade .tabset-pills}

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
 <thead>
  <tr>
   <th style="text-align:left;"> comparison </th>
   <th style="text-align:right;"> up </th>
   <th style="text-align:right;"> up % </th>
   <th style="text-align:right;"> down </th>
   <th style="text-align:right;"> down % </th>
   <th style="text-align:right;"> mean count cutoff </th>
   <th style="text-align:right;"> low counts </th>
   <th style="text-align:right;"> low counts % </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> sle - control </td>
   <td style="text-align:right;"> 4207 </td>
   <td style="text-align:right;"> 21.26 </td>
   <td style="text-align:right;"> 4568 </td>
   <td style="text-align:right;"> 23.09 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
<tfoot>
<tr><td style="padding: 0; " colspan="100%"><span style="font-style: italic;">Note: </span></td></tr>
<tr><td style="padding: 0; " colspan="100%">
<sup></sup> 19785 with nonzero total read count <br>Outliers:  0</td></tr>
</tfoot>
</table>

\newpage

## Volcano plot of gene expression

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/volcano_plot_of_gene_expression-1.png)<!-- -->

\newpage

## 25 genes with largest negative log-fold change in experimental group {.tabset .tabset-fade .tabset-pills}


<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption><b>sle - control: downregulated in sle - control</b></caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> gene </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> log2FoldChange </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> t </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> pvalue<sup>*</sup> </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> padj<sup>†</sup> </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> B </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">H3P6</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 100.00%">4.5</span> </td>
   <td style="text-align:right;"> -11.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">3.3e-24</span> </td>
   <td style="text-align:right;"> 48.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">TMEM189-UBE2V1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 97.78%">4.4</span> </td>
   <td style="text-align:right;"> -44.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">6e-143</span> </td>
   <td style="text-align:right;"> 290.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RPL21P16</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 84.44%">3.8</span> </td>
   <td style="text-align:right;"> -11.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.1e-23</span> </td>
   <td style="text-align:right;"> 47.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RPL22P1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 60.00%">2.7</span> </td>
   <td style="text-align:right;"> -8.2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">3.6e-14</span> </td>
   <td style="text-align:right;"> 24.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">TECRP1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 57.78%">2.6</span> </td>
   <td style="text-align:right;"> -10.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">4.5e-21</span> </td>
   <td style="text-align:right;"> 41.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">LYPD2</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 57.78%">2.6</span> </td>
   <td style="text-align:right;"> -15.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.6e-38</span> </td>
   <td style="text-align:right;"> 83.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RN7SL4P</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 57.78%">2.6</span> </td>
   <td style="text-align:right;"> -9.5 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">7.2e-18</span> </td>
   <td style="text-align:right;"> 33.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">SNORD10</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 55.56%">2.5</span> </td>
   <td style="text-align:right;"> -9.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.7e-16</span> </td>
   <td style="text-align:right;"> 30.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">PWAR5</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 55.56%">2.5</span> </td>
   <td style="text-align:right;"> -11.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">4.7e-24</span> </td>
   <td style="text-align:right;"> 47.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">KCNG1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 53.33%">2.4</span> </td>
   <td style="text-align:right;"> -15.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">7.3e-37</span> </td>
   <td style="text-align:right;"> 79.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">SNORA57</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 51.11%">2.3</span> </td>
   <td style="text-align:right;"> -8.1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.1e-13</span> </td>
   <td style="text-align:right;"> 23.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RPS28P7</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 51.11%">2.3</span> </td>
   <td style="text-align:right;"> -9.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2e-16</span> </td>
   <td style="text-align:right;"> 29.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RPL14P1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 48.89%">2.2</span> </td>
   <td style="text-align:right;"> -11.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">8.7e-23</span> </td>
   <td style="text-align:right;"> 45.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">SNORA48</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 48.89%">2.2</span> </td>
   <td style="text-align:right;"> -7.8 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">8.8e-13</span> </td>
   <td style="text-align:right;"> 21.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RPL10P9</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 46.67%">2.1</span> </td>
   <td style="text-align:right;"> -7.9 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">4.4e-13</span> </td>
   <td style="text-align:right;"> 21.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">SNORD33</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 46.67%">2.1</span> </td>
   <td style="text-align:right;"> -8.2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">5.6e-14</span> </td>
   <td style="text-align:right;"> 24.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">F8A3</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 46.67%">2.1</span> </td>
   <td style="text-align:right;"> -9.7 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2.2e-18</span> </td>
   <td style="text-align:right;"> 34.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">SNORD95</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 44.44%">2.0</span> </td>
   <td style="text-align:right;"> -18.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.3e-49</span> </td>
   <td style="text-align:right;"> 110.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">TUBB2A</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 40.00%">1.8</span> </td>
   <td style="text-align:right;"> -6.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(253, 231, 37, 1) !important;">1.9e-08</span> </td>
   <td style="text-align:right;"> 9.8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">SYN3</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 40.00%">1.8</span> </td>
   <td style="text-align:right;"> -11.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">5.1e-22</span> </td>
   <td style="text-align:right;"> 43.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">STEAP1B</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 40.00%">1.8</span> </td>
   <td style="text-align:right;"> -12.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">8.7e-28</span> </td>
   <td style="text-align:right;"> 57.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">TCL1A</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 40.00%">1.8</span> </td>
   <td style="text-align:right;"> -14.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.4e-32</span> </td>
   <td style="text-align:right;"> 69.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RNU1-1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 37.78%">1.7</span> </td>
   <td style="text-align:right;"> -11.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">3.9e-21</span> </td>
   <td style="text-align:right;"> 41.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RNU2-2P</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 37.78%">1.7</span> </td>
   <td style="text-align:right;"> -8.9 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">3.7e-16</span> </td>
   <td style="text-align:right;"> 29.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">IGHD</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 37.78%">1.7</span> </td>
   <td style="text-align:right;"> -13.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">9.1e-31</span> </td>
   <td style="text-align:right;"> 64.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">OR2T11</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 37.78%">1.7</span> </td>
   <td style="text-align:right;"> -24.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.3e-71</span> </td>
   <td style="text-align:right;"> 150.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">TMEM238</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: tomato; width: 37.78%">1.7</span> </td>
   <td style="text-align:right;"> -14.0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2.6e-34</span> </td>
   <td style="text-align:right;"> 73.0 </td>
  </tr>
</tbody>
<tfoot>
<tr><td style="padding: 0; " colspan="100%">
<sup>*</sup> Wald test p-values</td></tr>
<tr><td style="padding: 0; " colspan="100%">
<sup>†</sup> Benjamini–Hochberg adjusted value</td></tr>
</tfoot>
</table>

## 25 genes with largest positive log-fold change in experimental group {.tabset .tabset-fade .tabset-pills}


<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption><b>sle - control: upregulated in sle - control</b></caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> gene </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> log2FoldChange </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> F </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> pvalue<sup>*</sup> </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> padj<sup>†</sup> </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> logCPM </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">OTOF</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 100.00%">3.8</span> </td>
   <td style="text-align:right;"> 130 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">5.6e-23</span> </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">IFI27</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 92.11%">3.5</span> </td>
   <td style="text-align:right;"> 92 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.5e-17</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">IFI44L</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 78.95%">3.0</span> </td>
   <td style="text-align:right;"> 130 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">6.4e-24</span> </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RSAD2</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 71.05%">2.7</span> </td>
   <td style="text-align:right;"> 130 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">6e-23</span> </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">H2AC18</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 65.79%">2.5</span> </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(70, 7, 90, 1) !important;">4.4e-07</span> </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">IFI44</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 65.79%">2.5</span> </td>
   <td style="text-align:right;"> 160 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">4.2e-27</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">SIGLEC1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 65.79%">2.5</span> </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">8.3e-21</span> </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">USP18</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 60.53%">2.3</span> </td>
   <td style="text-align:right;"> 140 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2.2e-25</span> </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">ADAMTS2</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 57.89%">2.2</span> </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">6.9e-14</span> </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">CMPK2</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 57.89%">2.2</span> </td>
   <td style="text-align:right;"> 130 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">4.7e-24</span> </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">DEFA1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 57.89%">2.2</span> </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 4.0e-06 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(205, 225, 29, 1) !important;">2.5e-05</span> </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">IFIT1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 57.89%">2.2</span> </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">7.4e-19</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">ISG15</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 57.89%">2.2</span> </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">6.6e-20</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">OAS3</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 55.26%">2.1</span> </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">7.1e-21</span> </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">NRIR</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 52.63%">2.0</span> </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.8e-26</span> </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">EPSTI1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 50.00%">1.9</span> </td>
   <td style="text-align:right;"> 120 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">6.7e-22</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">IFIT3</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 50.00%">1.9</span> </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.4e-20</span> </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">LY6E</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 50.00%">1.9</span> </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.7e-19</span> </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">OAS1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 50.00%">1.9</span> </td>
   <td style="text-align:right;"> 130 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">5.2e-24</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">SPATS2L</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 50.00%">1.9</span> </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2.3e-26</span> </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">HERC5</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 47.37%">1.8</span> </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2.5e-19</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">MX1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 47.37%">1.8</span> </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">4.3e-19</span> </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">IFI6</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 44.74%">1.7</span> </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">6.2e-19</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">OASL</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 44.74%">1.7</span> </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">3.4e-20</span> </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">DDX60</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2.4e-26</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">ETV7</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">4e-19</span> </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">FBXO39</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2.8e-19</span> </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">LINC00487</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 90 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">2.6e-17</span> </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">MDGA1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.5e-08</span> </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">OAS2</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 130 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">9.2e-23</span> </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">PLSCR1</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 200 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">1.8e-32</span> </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RAP1GAP</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 4.3e-06 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(253, 231, 37, 1) !important;">2.7e-05</span> </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RN7SL605P</span> </td>
   <td style="text-align:left;color: black !important;"> <span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: lightgreen; width: 42.11%">1.6</span> </td>
   <td style="text-align:right;"> 240 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:left;"> <span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(68, 1, 84, 1) !important;">3.3e-38</span> </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
</tbody>
<tfoot>
<tr><td style="padding: 0; " colspan="100%">
<sup>*</sup> Wald test p-values</td></tr>
<tr><td style="padding: 0; " colspan="100%">
<sup>†</sup> Benjamini–Hochberg adjusted value</td></tr>
</tfoot>
</table>

\newpage


<!-- TODO: This is kind of dumb for the moment.  It only really does one group comparison, that is we only calculate -->
<!-- say the comparisons for each experimental group within Disease Class, not within Disease Class and within cluster and within race -->
<!-- but a lot more will have to change before we can do anything with that here -->
<!-- TODO: Add GSEA of top DEGs -->

```
## 
## 
## ## Expression of top DE genes for Disease_Class: sle - control {.tabset .tabset-fade .tabset-pills}
## 
## ### Sorted by Disease_Class 
## 
## ![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-13-1.png)<!-- -->
## 
## ### Sorted by cluster 
## 
## ![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-14-1.png)<!-- -->
## 
## ### Sorted by hierarchical_clustering 
## 
## ![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-15-1.png)<!-- -->
```

\newpage



## Pathway enrichment upgregulated genes for Disease_Class: sle - control {.tabset .tabset-fade .tabset-pills}

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

\newpage


## Pathway enrichment downgregulated genes for Disease_Class: sle - control {.tabset .tabset-fade .tabset-pills}

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

\newpage




## Expression of selected interferon-stimulated genes {.tabset .tabset-fade .tabset-pills}



### Sorted by Disease_Class


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-20-1.png)<!-- -->


### Sorted by cluster


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-21-1.png)<!-- -->


### Sorted by hierarchical clustering


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-22-1.png)<!-- -->


## Expression of individual genes in the M1.2, M3.4, and M5.12 gene modules {.tabset .tabset-fade .tabset-pills}



### Sorted by Disease_Class


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-23-1.png)<!-- -->


### Sorted by cluster


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-24-1.png)<!-- -->


### Sorted by hierarchical clustering


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

\newpage

# Sample Clustering

The eigenvalues for the first 100 principal components were used to calculate a random forest-based adjacency matrix, which was in turn used to calculate the gap statistic to determine the optimal *k*-clusters.


```{=html}
<div id="htmlwidget-173c4cc47a4edf65817a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-173c4cc47a4edf65817a">{"x":{"filter":"none","vertical":false,"data":[["1","2"],["control","sle"],[46,84],[6,68],[35,100],[null,27]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Disease_Class<\/th>\n      <th>2<\/th>\n      <th>3<\/th>\n      <th>4<\/th>\n      <th>1<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```
![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/Sample_Clustering_chart-2.png)<!-- -->

\newpage

# Clustering and variation of samples in reduced dimensional space {.tabset .tabset-fade .tabset-pills}




## PCA {.tabset .tabset-fade .tabset-pills}


### pca by Disease_Class 

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

### pca by ethnicity 

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

### pca by K-means cluster 

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

### pca by plate 

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

## UMAP {.tabset .tabset-fade .tabset-pills}


### umap by Disease_Class 

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

### umap by ethnicity 

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

### umap by K-means cluster 

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

### umap by plate 

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

\newpage



# Module score heatmaps {.tabset .tabset-fade .tabset-pills}


## Banchereau interferon modules scores {.tabset .tabset-fade .tabset-pills}


### Sorted by Disease_Class

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-34-1.png)<!-- -->

### Sorted by cluster

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

### Sorted by hierarchical clustering

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-36-1.png)<!-- -->

## Banchereau inflammation modules scores {.tabset .tabset-fade .tabset-pills}


### Sorted by Disease_Class

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-37-1.png)<!-- -->

### Sorted by cluster

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-38-1.png)<!-- -->

### Sorted by hierarchical clustering

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-39-1.png)<!-- -->

## Low-density granulocyte modules scores {.tabset .tabset-fade .tabset-pills}


### Sorted by Disease_Class

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-40-1.png)<!-- -->

### Sorted by cluster

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

### Sorted by hierarchical clustering

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-42-1.png)<!-- -->




## Scores for annotated Banchereau modules {.tabset .tabset-fade .tabset-pills}



### Sorted by Disease_Class


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-43-1.png)<!-- -->


### Sorted by cluster


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-44-1.png)<!-- -->


### Sorted by hierarchical clustering


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-45-1.png)<!-- -->


## Scores for all Banchereau modules {.tabset .tabset-fade .tabset-pills}



### Sorted by Disease_Class


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-46-1.png)<!-- -->


### Sorted by cluster


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-47-1.png)<!-- -->


### Sorted by hierarchical clustering


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

## Expression of the annotated Banchereau modules {.tabset .tabset-fade .tabset-pills}




### By Disease_Class {.tabset .tabset-fade .tabset-pills}



![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-49-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-50-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-51-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-52-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-53-1.png)<!-- -->


### By cluster {.tabset .tabset-fade .tabset-pills}



![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-54-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-55-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-56-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-57-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-58-1.png)<!-- -->

\newpage

# WGCNA



## Expression similarity dendrogram

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/WGCNA_Expression_similarity_dendrogram-1.png)<!-- -->

## Expression of identified modules {.tabset .tabset-fade .tabset-pills}



\newpage


### By Disease_Class {.tabset .tabset-fade .tabset-pills}



![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-59-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-60-1.png)<!-- -->


### By cluster {.tabset .tabset-fade .tabset-pills}



![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-61-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-62-1.png)<!-- -->

\newpage

## WGCNA eigenvalues {.tabset .tabset-fade .tabset-pills}



\newpage


### Sorted by Disease_Class


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-63-1.png)<!-- -->


### Sorted by cluster


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-64-1.png)<!-- -->


### Sorted by hierarchical clustering


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-65-1.png)<!-- -->


\newpage

## Gene ontology GSEA of WGCNA modules




![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-66-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-67-1.png)<!-- -->

![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/unnamed-chunk-68-1.png)<!-- -->

\newpage

# Use of gene expression modules in classifying samples.





## WGCNA scores {.tabset .tabset-fade .tabset-pills}


### Sorted by Disease_Class


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/wgcna_gsea_plot_WGCNA_Disease_Class-1.png)<!-- -->

### Sorted by cluster


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/wgcna_gsea_plot_WGCNA_cluster-1.png)<!-- -->


## Banchereau scores {.tabset .tabset-fade .tabset-pills}


### Sorted by Disease_Class


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/wgcna_gsea_plot_Banchereau_Disease_Class-1.png)<!-- -->

### Sorted by cluster


![](/home/rstudio/workspace/rnaseq_targets_pipeline/reports/report_files/figure-html/wgcna_gsea_plot_Banchereau_cluster-1.png)<!-- -->


```
## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 4.1.1 (2021-08-10)
##  os       Ubuntu 20.04.3 LTS          
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  ctype    en_US.UTF-8                 
##  tz       Etc/UTC                     
##  date     2021-09-29                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  package          * version    date       lib source        
##  abind              1.4-5      2016-07-21 [1] RSPM (R 4.1.1)
##  annotate           1.70.0     2021-05-19 [1] RSPM (R 4.1.1)
##  AnnotationDbi      1.54.1     2021-06-08 [1] RSPM (R 4.1.1)
##  ape                5.5        2021-04-25 [1] RSPM (R 4.1.1)
##  aplot              0.1.0      2021-09-03 [1] RSPM (R 4.1.1)
##  ash                1.0-15     2015-09-01 [1] RSPM (R 4.1.0)
##  assertthat         0.2.1      2019-03-21 [1] RSPM (R 4.1.1)
##  backports          1.2.1      2020-12-09 [1] RSPM (R 4.1.1)
##  base64enc          0.1-3      2015-07-28 [1] RSPM (R 4.1.1)
##  beeswarm           0.4.0      2021-06-01 [1] RSPM (R 4.1.1)
##  Biobase            2.52.0     2021-05-19 [1] RSPM (R 4.1.1)
##  BiocGenerics       0.38.0     2021-05-19 [1] RSPM (R 4.1.1)
##  BiocParallel       1.26.2     2021-08-22 [1] RSPM (R 4.1.1)
##  Biostrings         2.60.2     2021-08-05 [1] RSPM (R 4.1.1)
##  bit                4.0.4      2020-08-04 [1] RSPM (R 4.1.0)
##  bit64              4.0.5      2020-08-30 [1] RSPM (R 4.1.0)
##  bitops             1.0-7      2021-04-24 [1] RSPM (R 4.1.1)
##  blob               1.2.2      2021-07-23 [1] RSPM (R 4.1.1)
##  broom              0.7.9      2021-07-27 [1] RSPM (R 4.1.1)
##  bslib              0.3.0      2021-09-02 [1] RSPM (R 4.1.1)
##  cachem             1.0.6      2021-08-19 [1] RSPM (R 4.1.0)
##  Cairo              1.5-12.2   2020-07-07 [1] RSPM (R 4.1.1)
##  callr              3.7.0      2021-04-20 [1] RSPM (R 4.1.0)
##  car                3.0-11     2021-06-27 [1] RSPM (R 4.1.1)
##  carData            3.0-4      2020-05-22 [1] RSPM (R 4.1.1)
##  cellranger         1.1.0      2016-07-27 [1] RSPM (R 4.1.1)
##  checkmate          2.0.0      2020-02-06 [1] RSPM (R 4.1.1)
##  circlize           0.4.13     2021-06-09 [1] RSPM (R 4.1.0)
##  class              7.3-19     2021-05-03 [2] CRAN (R 4.1.1)
##  cli                3.0.1      2021-07-17 [1] RSPM (R 4.1.0)
##  clue               0.3-59     2021-04-16 [1] RSPM (R 4.1.0)
##  cluster            2.1.2      2021-04-17 [2] CRAN (R 4.1.1)
##  codetools          0.2-18     2020-11-04 [2] CRAN (R 4.1.1)
##  colorspace         2.0-2      2021-06-24 [1] RSPM (R 4.1.1)
##  ComplexHeatmap     2.8.0      2021-05-19 [1] Bioconductor  
##  cowplot            1.1.1      2020-12-30 [1] RSPM (R 4.1.1)
##  crayon             1.4.1      2021-02-08 [1] RSPM (R 4.1.0)
##  crosstalk          1.1.1      2021-01-12 [1] RSPM (R 4.1.1)
##  curl               4.3.2      2021-06-23 [1] RSPM (R 4.1.0)
##  data.table         1.14.0     2021-02-21 [1] RSPM (R 4.1.1)
##  DBI                1.1.1      2021-01-15 [1] RSPM (R 4.1.1)
##  dials              0.0.10     2021-09-10 [1] RSPM (R 4.1.0)
##  DiceDesign         1.9        2021-02-13 [1] RSPM (R 4.1.0)
##  digest             0.6.27     2020-10-24 [1] RSPM (R 4.1.0)
##  DO.db              2.9        2021-09-17 [1] RSPM (R 4.1.1)
##  doParallel         1.0.16     2020-10-16 [1] RSPM (R 4.1.1)
##  DOSE               3.18.2     2021-08-17 [1] RSPM (R 4.1.1)
##  dplyr            * 1.0.7      2021-06-18 [1] RSPM (R 4.1.1)
##  DT                 0.19       2021-09-02 [1] RSPM (R 4.1.1)
##  dynamicTreeCut     1.63-1     2016-03-11 [1] RSPM (R 4.1.1)
##  ellipsis           0.3.2      2021-04-29 [1] RSPM (R 4.1.0)
##  EnhancedVolcano    1.10.0     2021-05-19 [1] Bioconductor  
##  enrichplot         1.12.2     2021-07-01 [1] RSPM (R 4.1.1)
##  evaluate           0.14       2019-05-28 [1] RSPM (R 4.1.0)
##  extrafont          0.17       2014-12-08 [1] RSPM (R 4.1.0)
##  extrafontdb        1.0        2012-06-11 [1] RSPM (R 4.1.0)
##  fansi              0.5.0      2021-05-25 [1] RSPM (R 4.1.0)
##  farver             2.1.0      2021-02-28 [1] RSPM (R 4.1.1)
##  fastcluster        1.2.3      2021-05-24 [1] RSPM (R 4.1.1)
##  fastmap            1.1.0      2021-01-25 [1] RSPM (R 4.1.0)
##  fastmatch          1.1-3      2021-07-23 [1] RSPM (R 4.1.1)
##  fgsea              1.18.0     2021-05-19 [1] RSPM (R 4.1.1)
##  flextable        * 0.6.8      2021-09-06 [1] RSPM (R 4.1.0)
##  forcats            0.5.1      2021-01-27 [1] RSPM (R 4.1.1)
##  foreach            1.5.1      2020-10-15 [1] RSPM (R 4.1.1)
##  foreign            0.8-81     2020-12-22 [2] CRAN (R 4.1.1)
##  formattable      * 0.2.1      2021-01-07 [1] RSPM (R 4.1.1)
##  Formula            1.2-4      2020-10-16 [1] RSPM (R 4.1.1)
##  fs                 1.5.0      2020-07-31 [1] RSPM (R 4.1.0)
##  furrr              0.2.3      2021-06-25 [1] RSPM (R 4.1.1)
##  future             1.22.1     2021-08-25 [1] RSPM (R 4.1.1)
##  future.apply       1.8.1      2021-08-10 [1] RSPM (R 4.1.0)
##  gdtools            0.2.3      2021-01-06 [1] RSPM (R 4.1.0)
##  genefilter       * 1.74.0     2021-05-19 [1] RSPM (R 4.1.1)
##  generics           0.1.0      2020-10-31 [1] RSPM (R 4.1.1)
##  GenomeInfoDb       1.28.4     2021-09-05 [1] RSPM (R 4.1.1)
##  GenomeInfoDbData   1.2.6      2021-09-17 [1] RSPM (R 4.1.1)
##  GetoptLong         1.0.5      2020-12-15 [1] RSPM (R 4.1.0)
##  ggalt              0.4.0      2017-02-15 [1] RSPM (R 4.1.0)
##  ggbeeswarm         0.6.0      2017-08-07 [1] RSPM (R 4.1.1)
##  ggforce            0.3.3      2021-03-05 [1] RSPM (R 4.1.1)
##  ggfun              0.0.4      2021-09-17 [1] RSPM (R 4.1.1)
##  ggnewscale         0.4.5      2021-01-11 [1] RSPM (R 4.1.0)
##  ggplot2          * 3.3.5      2021-06-25 [1] RSPM (R 4.1.1)
##  ggplotify          0.1.0      2021-09-02 [1] RSPM (R 4.1.1)
##  ggpubr           * 0.4.0      2020-06-27 [1] RSPM (R 4.1.1)
##  ggraph             2.0.5      2021-02-23 [1] RSPM (R 4.1.1)
##  ggrastr            0.2.3      2021-03-01 [1] RSPM (R 4.1.1)
##  ggrepel            0.9.1      2021-01-15 [1] RSPM (R 4.1.1)
##  ggsignif           0.6.3      2021-09-09 [1] RSPM (R 4.1.1)
##  ggtree             3.0.4      2021-08-22 [1] RSPM (R 4.1.1)
##  GlobalOptions      0.1.2      2020-06-10 [1] RSPM (R 4.1.0)
##  globals            0.14.0     2020-11-22 [1] RSPM (R 4.1.1)
##  glue               1.4.2      2020-08-27 [1] RSPM (R 4.1.0)
##  GO.db              3.13.0     2021-09-17 [1] RSPM (R 4.1.1)
##  GOSemSim           2.18.1     2021-07-29 [1] RSPM (R 4.1.1)
##  gower              0.2.2      2020-06-23 [1] RSPM (R 4.1.0)
##  GPfit              1.0-8      2019-02-08 [1] RSPM (R 4.1.0)
##  graphlayouts       0.7.1      2020-10-26 [1] RSPM (R 4.1.1)
##  gridExtra          2.3        2017-09-09 [1] RSPM (R 4.1.1)
##  gridGraphics       0.5-1      2020-12-13 [1] RSPM (R 4.1.1)
##  gtable             0.3.0      2019-03-25 [1] RSPM (R 4.1.1)
##  hardhat            0.1.6      2021-07-14 [1] RSPM (R 4.1.0)
##  haven              2.4.3      2021-08-04 [1] RSPM (R 4.1.1)
##  here             * 1.0.1      2020-12-13 [1] RSPM (R 4.1.1)
##  highr              0.9        2021-04-16 [1] RSPM (R 4.1.0)
##  Hmisc              4.5-0      2021-02-28 [1] RSPM (R 4.1.1)
##  hms                1.1.0      2021-05-17 [1] RSPM (R 4.1.1)
##  htmlTable          2.2.1      2021-05-18 [1] RSPM (R 4.1.1)
##  htmltools          0.5.2      2021-08-25 [1] RSPM (R 4.1.1)
##  htmlwidgets        1.5.4      2021-09-08 [1] RSPM (R 4.1.1)
##  httr               1.4.2      2020-07-20 [1] RSPM (R 4.1.0)
##  igraph             1.2.6      2020-10-06 [1] RSPM (R 4.1.1)
##  impute             1.66.0     2021-05-19 [1] RSPM (R 4.1.1)
##  ipred              0.9-12     2021-09-15 [1] RSPM (R 4.1.1)
##  IRanges            2.26.0     2021-05-19 [1] RSPM (R 4.1.1)
##  iterators          1.0.13     2020-10-15 [1] RSPM (R 4.1.1)
##  janitor          * 2.1.0      2021-01-05 [1] RSPM (R 4.1.1)
##  jpeg               0.1-9      2021-07-24 [1] RSPM (R 4.1.1)
##  jquerylib          0.1.4      2021-04-26 [1] RSPM (R 4.1.1)
##  jsonlite           1.7.2      2020-12-09 [1] RSPM (R 4.1.0)
##  kableExtra       * 1.3.4      2021-02-20 [1] RSPM (R 4.1.1)
##  KEGGREST           1.32.0     2021-05-19 [1] RSPM (R 4.1.1)
##  KernSmooth         2.23-20    2021-05-03 [2] CRAN (R 4.1.1)
##  knitr            * 1.34       2021-09-09 [1] RSPM (R 4.1.0)
##  labeling           0.4.2      2020-10-20 [1] RSPM (R 4.1.1)
##  lattice            0.20-44    2021-05-02 [2] CRAN (R 4.1.1)
##  latticeExtra       0.6-29     2019-12-19 [1] RSPM (R 4.1.1)
##  lava               1.6.10     2021-09-02 [1] RSPM (R 4.1.0)
##  lazyeval           0.2.2      2019-03-15 [1] RSPM (R 4.1.1)
##  lhs                1.1.3      2021-09-08 [1] RSPM (R 4.1.0)
##  lifecycle          1.0.0      2021-02-15 [1] RSPM (R 4.1.0)
##  listenv            0.8.0      2019-12-05 [1] RSPM (R 4.1.1)
##  lubridate          1.7.10     2021-02-26 [1] RSPM (R 4.1.1)
##  magrittr         * 2.0.1      2020-11-17 [1] RSPM (R 4.1.0)
##  maps               3.4.0      2021-09-25 [1] RSPM (R 4.1.0)
##  MASS               7.3-54     2021-05-03 [2] CRAN (R 4.1.1)
##  Matrix             1.3-4      2021-06-01 [2] CRAN (R 4.1.1)
##  matrixStats        0.61.0     2021-09-17 [1] RSPM (R 4.1.1)
##  memoise            2.0.0      2021-01-26 [1] RSPM (R 4.1.0)
##  munsell            0.5.0      2018-06-12 [1] RSPM (R 4.1.1)
##  nlme               3.1-153    2021-09-07 [2] RSPM (R 4.1.0)
##  nnet               7.3-16     2021-05-03 [2] CRAN (R 4.1.1)
##  officer            0.4.0      2021-09-06 [1] RSPM (R 4.1.0)
##  openxlsx           4.2.4      2021-06-16 [1] RSPM (R 4.1.1)
##  paletteer        * 1.4.0      2021-07-20 [1] RSPM (R 4.1.1)
##  parallelly         1.28.1     2021-09-09 [1] RSPM (R 4.1.1)
##  parsnip            0.1.7      2021-07-21 [1] RSPM (R 4.1.0)
##  patchwork          1.1.1      2020-12-17 [1] RSPM (R 4.1.1)
##  pillar             1.6.2      2021-07-29 [1] RSPM (R 4.1.0)
##  pkgconfig          2.0.3      2019-09-22 [1] RSPM (R 4.1.0)
##  plyr               1.8.6      2020-03-03 [1] RSPM (R 4.1.1)
##  png                0.1-7      2013-12-03 [1] RSPM (R 4.1.1)
##  polyclip           1.10-0     2019-03-14 [1] RSPM (R 4.1.1)
##  preprocessCore     1.54.0     2021-05-19 [1] RSPM (R 4.1.1)
##  prismatic          1.0.0      2021-01-05 [1] RSPM (R 4.1.1)
##  pROC               1.18.0     2021-09-03 [1] RSPM (R 4.1.0)
##  processx           3.5.2      2021-04-30 [1] RSPM (R 4.1.0)
##  prodlim            2019.11.13 2019-11-17 [1] RSPM (R 4.1.0)
##  proj4              1.0-10.1   2021-01-26 [1] RSPM (R 4.1.0)
##  ps                 1.6.0      2021-02-28 [1] RSPM (R 4.1.0)
##  purrr            * 0.3.4      2020-04-17 [1] RSPM (R 4.1.0)
##  qvalue             2.24.0     2021-05-19 [1] RSPM (R 4.1.1)
##  R6                 2.5.1      2021-08-19 [1] RSPM (R 4.1.0)
##  ranger             0.13.1     2021-07-14 [1] RSPM (R 4.1.0)
##  RColorBrewer       1.1-2      2014-12-07 [1] RSPM (R 4.1.1)
##  Rcpp               1.0.7      2021-07-07 [1] RSPM (R 4.1.1)
##  RCurl              1.98-1.5   2021-09-17 [1] RSPM (R 4.1.1)
##  readxl             1.3.1      2019-03-13 [1] RSPM (R 4.1.1)
##  recipes            0.1.16     2021-04-16 [1] RSPM (R 4.1.0)
##  rematch2           2.1.2      2020-05-01 [1] RSPM (R 4.1.0)
##  reshape2           1.4.4      2020-04-09 [1] RSPM (R 4.1.1)
##  rio                0.5.27     2021-06-21 [1] RSPM (R 4.1.1)
##  rjson              0.2.20     2018-06-08 [1] RSPM (R 4.1.1)
##  rlang            * 0.4.11     2021-04-30 [1] RSPM (R 4.1.0)
##  rmarkdown        * 2.11       2021-09-14 [1] RSPM (R 4.1.1)
##  rpart              4.1-15     2019-04-12 [2] CRAN (R 4.1.1)
##  rprojroot          2.0.2      2020-11-15 [1] RSPM (R 4.1.0)
##  rsample            0.1.0      2021-05-08 [1] RSPM (R 4.1.0)
##  RSQLite            2.2.8      2021-08-21 [1] RSPM (R 4.1.1)
##  rstatix            0.7.0      2021-02-13 [1] RSPM (R 4.1.1)
##  rstudioapi         0.13       2020-11-12 [1] RSPM (R 4.1.0)
##  Rttf2pt1           1.3.9      2021-07-22 [1] RSPM (R 4.1.0)
##  rvest              1.0.1      2021-07-26 [1] RSPM (R 4.1.1)
##  S4Vectors          0.30.0     2021-05-19 [1] RSPM (R 4.1.1)
##  sass               0.4.0      2021-05-12 [1] RSPM (R 4.1.1)
##  scales             1.1.1      2020-05-11 [1] RSPM (R 4.1.1)
##  scatterpie         0.1.7      2021-08-20 [1] RSPM (R 4.1.1)
##  sessioninfo        1.1.1      2018-11-05 [1] RSPM (R 4.1.0)
##  shadowtext         0.0.9      2021-09-19 [1] RSPM (R 4.1.1)
##  shape              1.4.6      2021-05-19 [1] RSPM (R 4.1.0)
##  snakecase          0.11.0     2019-05-25 [1] RSPM (R 4.1.1)
##  stringi            1.7.4      2021-08-25 [1] RSPM (R 4.1.0)
##  stringr          * 1.4.0      2019-02-10 [1] RSPM (R 4.1.0)
##  survival           3.2-13     2021-08-24 [2] RSPM (R 4.1.0)
##  svglite            2.0.0      2021-02-20 [1] RSPM (R 4.1.1)
##  systemfonts        1.0.2      2021-05-11 [1] RSPM (R 4.1.1)
##  tarchetypes      * 0.3.0      2021-08-04 [1] RSPM (R 4.1.0)
##  targets          * 0.7.0      2021-08-19 [1] RSPM (R 4.1.1)
##  tibble           * 3.1.4      2021-08-25 [1] RSPM (R 4.1.0)
##  tidygraph          1.2.0      2020-05-12 [1] RSPM (R 4.1.1)
##  tidyr            * 1.1.3      2021-03-03 [1] RSPM (R 4.1.1)
##  tidyselect       * 1.1.1      2021-04-30 [1] RSPM (R 4.1.1)
##  tidytree           0.3.5      2021-09-08 [1] RSPM (R 4.1.1)
##  timeDate           3043.102   2018-02-21 [1] RSPM (R 4.1.0)
##  treeio             1.16.2     2021-08-17 [1] RSPM (R 4.1.1)
##  tune               0.1.6      2021-07-21 [1] RSPM (R 4.1.0)
##  tweenr             1.0.2      2021-03-23 [1] RSPM (R 4.1.1)
##  utf8               1.2.2      2021-07-24 [1] RSPM (R 4.1.0)
##  uuid               0.1-4      2020-02-26 [1] RSPM (R 4.1.1)
##  vctrs              0.3.8      2021-04-29 [1] RSPM (R 4.1.0)
##  vip                0.3.2      2020-12-17 [1] RSPM (R 4.1.0)
##  vipor              0.4.5      2017-03-22 [1] RSPM (R 4.1.1)
##  viridis            0.6.1      2021-05-11 [1] RSPM (R 4.1.1)
##  viridisLite        0.4.0      2021-04-13 [1] RSPM (R 4.1.1)
##  webshot            0.5.2      2019-11-22 [1] RSPM (R 4.1.1)
##  WGCNA              1.70-3     2021-02-28 [1] RSPM (R 4.1.1)
##  withr              2.4.2      2021-04-18 [1] RSPM (R 4.1.0)
##  workflows          0.2.3      2021-07-16 [1] RSPM (R 4.1.0)
##  xfun               0.26       2021-09-14 [1] RSPM (R 4.1.1)
##  XML                3.99-0.8   2021-09-17 [1] RSPM (R 4.1.1)
##  xml2               1.3.2      2020-04-23 [1] RSPM (R 4.1.0)
##  xtable             1.8-4      2019-04-21 [1] RSPM (R 4.1.1)
##  XVector            0.32.0     2021-05-19 [1] RSPM (R 4.1.1)
##  yaml               2.2.1      2020-02-01 [1] RSPM (R 4.1.0)
##  yardstick          0.0.8      2021-03-28 [1] RSPM (R 4.1.0)
##  yulab.utils        0.0.2      2021-08-16 [1] RSPM (R 4.1.1)
##  zip                2.2.0      2021-05-31 [1] RSPM (R 4.1.0)
##  zlibbioc           1.38.0     2021-05-19 [1] RSPM (R 4.1.1)
## 
## [1] /usr/local/lib/R/site-library
## [2] /usr/local/lib/R/library
```
