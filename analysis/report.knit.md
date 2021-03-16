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
title: "Initial COVID PCV samples RNAseq Analysis"
author: "Miles Smith"
date: "15 March, 2021"
output: 
  html_document:
    highlight: "pygments"
    toc: TRUE
    toc_depth: 3
    toc_float: TRUE
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





```{=html}
<template id="264ba4c8-fb75-486f-b6a9-ca112b9b0f07"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 5px;
  margin-bottom: 5px;
  table-layout: fixed;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-2b9ea57c{table-layout:auto;border-collapse:collapse;width:100%;}.cl-2b9639b4{font-family:'DejaVu Sans';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-2b96618c{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-2b9661c8{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-2b96b146{width:133pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2b96b178{width:57pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2b96b18c{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2b96b1a0{width:133pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2b96b1aa{width:57pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2b96b1be{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2b96b1d2{width:133pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2b96b1e6{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2b96b1f0{width:57pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-2b9ea57c'><thead><tr style="overflow-wrap:break-word;"><td class="cl-2b96b1d2"><p class="cl-2b96618c"><span class="cl-2b9639b4">Disease classification</span></p></td><td class="cl-2b96b1e6"><p class="cl-2b9661c8"><span class="cl-2b9639b4">Number</span></p></td><td class="cl-2b96b1f0"><p class="cl-2b9661c8"><span class="cl-2b9639b4">Percent</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-2b96b146"><p class="cl-2b96618c"><span class="cl-2b9639b4">control</span></p></td><td class="cl-2b96b18c"><p class="cl-2b9661c8"><span class="cl-2b9639b4">87</span></p></td><td class="cl-2b96b178"><p class="cl-2b9661c8"><span class="cl-2b9639b4">48.1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2b96b1a0"><p class="cl-2b96618c"><span class="cl-2b9639b4">infected</span></p></td><td class="cl-2b96b1be"><p class="cl-2b9661c8"><span class="cl-2b9639b4">94</span></p></td><td class="cl-2b96b1aa"><p class="cl-2b9661c8"><span class="cl-2b9639b4">51.9</span></p></td></tr></tbody></table></div></template>
<div id="0156177c-d79d-47db-a5e3-c0f17f99f7db"></div>
<script>
var dest = document.getElementById("0156177c-d79d-47db-a5e3-c0f17f99f7db");
var template = document.getElementById("264ba4c8-fb75-486f-b6a9-ca112b9b0f07");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;"
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
fantome.appendChild(templateContent);
</script>

```

\newpage

# Differential Gene Expression {.tabset .tabset-fade .tabset-pills}
```{=html}
<template id="a6b0061e-59a3-4406-946b-3dd2479fba62"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 5px;
  margin-bottom: 5px;
  table-layout: fixed;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-2bce88d2{table-layout:auto;border-collapse:collapse;width:100%;}.cl-2bc6cb4c{font-family:'DejaVu Sans';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-2bc6cb92{font-family:'DejaVu Sans';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-2bc6e80c{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-2bc6e83e{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-2bc72b0a{width:119pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72b46{width:93pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72b64{width:103pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72b78{width:120pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72b8c{width:94pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72ba0{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72baa{width:118pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72bbe{width:118pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72bc8{width:119pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72bdc{width:93pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72bf0{width:103pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72c04{width:120pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72c0e{width:94pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72c2c{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72c36{width:118pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2bc72c40{width:118pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-2bce88d2'><thead><tr style="overflow-wrap:break-word;"><td class="cl-2bc72c36"><p class="cl-2bc6e80c"><span class="cl-2bc6cb4c">Comparison</span></p></td><td class="cl-2bc72bf0"><p class="cl-2bc6e83e"><span class="cl-2bc6cb4c"># Up<br>regulated</span></p></td><td class="cl-2bc72c2c"><p class="cl-2bc6e83e"><span class="cl-2bc6cb4c">% Up<br>regulated</span></p></td><td class="cl-2bc72bc8"><p class="cl-2bc6e83e"><span class="cl-2bc6cb4c"># Down<br>regulated</span></p></td><td class="cl-2bc72c04"><p class="cl-2bc6e83e"><span class="cl-2bc6cb4c">% Down<br>regulated</span></p></td><td class="cl-2bc72c40"><p class="cl-2bc6e83e"><span class="cl-2bc6cb4c">Mean count<br>cutoff</span></p></td><td class="cl-2bc72bdc"><p class="cl-2bc6e83e"><span class="cl-2bc6cb4c"># Low<br>counts</span></p></td><td class="cl-2bc72c0e"><p class="cl-2bc6e83e"><span class="cl-2bc6cb4c">% Low<br>counts</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-2bc72baa"><p class="cl-2bc6e80c"><span class="cl-2bc6cb92">infected vs control</span></p></td><td class="cl-2bc72b64"><p class="cl-2bc6e83e"><span class="cl-2bc6cb92">4,207</span></p></td><td class="cl-2bc72ba0"><p class="cl-2bc6e83e"><span class="cl-2bc6cb92">9.96</span></p></td><td class="cl-2bc72b0a"><p class="cl-2bc6e83e"><span class="cl-2bc6cb92">3,534</span></p></td><td class="cl-2bc72b78"><p class="cl-2bc6e83e"><span class="cl-2bc6cb92">8.36</span></p></td><td class="cl-2bc72bbe"><p class="cl-2bc6e83e"><span class="cl-2bc6cb92">0</span></p></td><td class="cl-2bc72b46"><p class="cl-2bc6e83e"><span class="cl-2bc6cb92">12,289</span></p></td><td class="cl-2bc72b8c"><p class="cl-2bc6e83e"><span class="cl-2bc6cb92">29.08</span></p></td></tr></tbody></table></div></template>
<div id="9df56bd0-1b81-474d-b40a-ff852ba66208"></div>
<script>
var dest = document.getElementById("9df56bd0-1b81-474d-b40a-ff852ba66208");
var template = document.getElementById("a6b0061e-59a3-4406-946b-3dd2479fba62");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;"
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
fantome.appendChild(templateContent);
</script>

```

\newpage

## 25 genes with largest negative log-fold change in experimental group {.tabset .tabset-fade .tabset-pills}


```{=html}
<template id="1e5329d5-47c9-4d3c-ac91-825a0c43baf7"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 5px;
  margin-bottom: 5px;
  table-layout: fixed;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-2c09fcfa{table-layout:auto;border-collapse:collapse;width:100%;}.cl-2bffc82a{font-family:'DejaVu Sans';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-2bffc870{font-family:'DejaVu Sans';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-2bffc884{font-family:'DejaVu Sans';font-size:7pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;top:3pt;}.cl-2bffc898{font-family:'DejaVu Sans';font-size:7pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;bottom:3pt;}.cl-2bffc8f2{font-family:'DejaVu Sans';font-size:7pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;bottom:3pt;}.cl-2bffe382{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-2bffe3b4{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-2c00660e{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c00665e{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006672{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006686{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c00669a{width:83pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0066ae{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0066c2{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0066d6{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0066ea{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006708{width:83pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c00671c{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c00673a{width:83pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c00674e{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006762{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006776{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c00678a{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006794{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0067a8{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0067b2{width:83pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0067c6{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0067da{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c0067e4{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006802{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006816{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c006820{width:83pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-2c09fcfa'><thead><tr style="overflow-wrap:break-word;"><td  colspan="5"class="cl-2c0067b2"><p class="cl-2bffe382"><span class="cl-2bffc82a">infected_vs_control</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c006820"><p class="cl-2bffe382"><span class="cl-2bffc870">Gene</span></p></td><td class="cl-2c006816"><p class="cl-2bffe3b4"><span class="cl-2bffc870">Log</span><span class="cl-2bffc884">2</span><span class="cl-2bffc870">-fold change</span></p></td><td class="cl-2c006802"><p class="cl-2bffe3b4"><span class="cl-2bffc870">Log-fold change<br>standard error</span></p></td><td class="cl-2c0067e4"><p class="cl-2bffe3b4"><span class="cl-2bffc870">P-value</span><span class="cl-2bffc898">*</span></p></td><td class="cl-2c0067da"><p class="cl-2bffe3b4"><span class="cl-2bffc870">Adjusted<br>P-value</span><span class="cl-2bffc898">†</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC116533.1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">13</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">3.60</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">6.7e-167</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.8e-164</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC090227.1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">12</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">3.20</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.6e-53</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.7e-51</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">SCARNA17</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">11</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">3.30</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.2e-10</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.4e-09</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AP000944.4</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">9.8</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">3.00</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.3e-09</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.0e-08</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC005258.1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">8.9</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.00</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">7.1e-35</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.9e-33</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">H3F3AP4</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">7.7</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.81</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">9.5e-32</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5.9e-30</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">RXYLT1-AS1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">7.6</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.80</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.1e-03</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">6.5e-03</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">MIR7161</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">7.5</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.80</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.9e-03</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.5e-02</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AL355916.3</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">7.4</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.75</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5.6e-60</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">6.7e-58</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC016573.1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">7</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.40</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.1e-11</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.6e-10</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">LINC00395</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">6.5</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.30</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">8.5e-15</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.5e-13</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC105383.1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">6.3</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.30</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">6.2e-131</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.6e-128</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC104211.4</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5.9</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.40</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">3.4e-05</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.9e-04</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">RPL22P1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5.6</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.59</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.2e-28</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">6.5e-27</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">PPIAP11</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5.5</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.26</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">9.4e-126</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.4e-123</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC093843.2</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5.4</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.40</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.1e-56</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.2e-54</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">LINC02016</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5.4</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.50</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.9e-03</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.5e-02</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">RN7SL586P</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5.4</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.58</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.2e-38</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.6e-36</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AP000873.3</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">5</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.64</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.7e-23</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">7.9e-22</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC093765.2</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.9</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.35</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.2e-72</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">3.1e-70</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">MTCO1P40</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.9</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.20</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.4e-08</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.2e-07</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC083923.1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.8</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.10</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.0e-04</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.6e-03</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC093772.1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.7</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.34</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">3.0e-67</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.0e-65</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">PPIAP79</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.7</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.27</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">2.1e-82</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">3.4e-80</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c00669a"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC007389.1</span></p></td><td class="cl-2c006686"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.4</span></p></td><td class="cl-2c00660e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.25</span></p></td><td class="cl-2c006672"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">6.3e-72</span></p></td><td class="cl-2c00665e"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">8.8e-70</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c006708"><p class="cl-2bffe382"><span class="cl-2bffc82a">AC124242.3</span></p></td><td class="cl-2c0066ae"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">4.4</span></p></td><td class="cl-2c0066ea"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">0.20</span></p></td><td class="cl-2c0066c2"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">8.7e-111</span></p></td><td class="cl-2c0066d6"><p class="cl-2bffe3b4"><span class="cl-2bffc82a">1.9e-108</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="5"class="cl-2c00673a"><p class="cl-2bffe382"><span class="cl-2bffc8f2">*</span><span class="cl-2bffc82a">Wald test p-values</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="5"class="cl-2c00673a"><p class="cl-2bffe382"><span class="cl-2bffc8f2">†</span><span class="cl-2bffc82a">Benjamini–Hochberg adjusted value</span></p></td></tr></tfoot></table></div></template>
<div id="bb3abb95-49b3-4863-9464-918d86aaf7b9"></div>
<script>
var dest = document.getElementById("bb3abb95-49b3-4863-9464-918d86aaf7b9");
var template = document.getElementById("1e5329d5-47c9-4d3c-ac91-825a0c43baf7");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;"
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
fantome.appendChild(templateContent);
</script>

```

## 25 genes with largest positive log-fold change in experimental group {.tabset .tabset-fade .tabset-pills}


```{=html}
<template id="0a79c07c-ba2b-4576-b32a-7f98b7256987"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 5px;
  margin-bottom: 5px;
  table-layout: fixed;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-2c489186{table-layout:auto;border-collapse:collapse;width:100%;}.cl-2c3edbe6{font-family:'DejaVu Sans';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-2c3edc2c{font-family:'DejaVu Sans';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-2c3edc36{font-family:'DejaVu Sans';font-size:7pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;top:3pt;}.cl-2c3edc4a{font-family:'DejaVu Sans';font-size:7pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;bottom:3pt;}.cl-2c3edc54{font-family:'DejaVu Sans';font-size:7pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;bottom:3pt;}.cl-2c3ef7f2{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-2c3ef824{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-2c3f6b6a{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6ba6{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6bba{width:90pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6bc4{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6bd8{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6bec{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c00{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c14{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c1e{width:90pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c32{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c50{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c64{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c6e{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c82{width:90pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6c8c{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6ca0{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6cb4{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6cc8{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6cdc{width:90pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6ce6{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6cfa{width:111pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6d04{width:190pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6d18{width:90pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6d2c{width:64pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-2c3f6d36{width:104pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-2c489186'><thead><tr style="overflow-wrap:break-word;"><td  colspan="5"class="cl-2c3f6cdc"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">infected_vs_control</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6d18"><p class="cl-2c3ef7f2"><span class="cl-2c3edc2c">Gene</span></p></td><td class="cl-2c3f6d36"><p class="cl-2c3ef824"><span class="cl-2c3edc2c">Log</span><span class="cl-2c3edc36">2</span><span class="cl-2c3edc2c">-fold change</span></p></td><td class="cl-2c3f6d04"><p class="cl-2c3ef824"><span class="cl-2c3edc2c">Log-fold change<br>standard error</span></p></td><td class="cl-2c3f6d2c"><p class="cl-2c3ef824"><span class="cl-2c3edc2c">P-value</span><span class="cl-2c3edc4a">*</span></p></td><td class="cl-2c3f6cfa"><p class="cl-2c3ef824"><span class="cl-2c3edc2c">Adjusted<br>P-value</span><span class="cl-2c3edc4a">†</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">CTBP2P3</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.8</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.70</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.3e-12</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">5.4e-11</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">RN7SL571P</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">7.6</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.70</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.2e-05</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.5e-04</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AC073071.1</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">7.2</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.50</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.7e-05</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.2e-04</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AC097110.1</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">5.4</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.40</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.7e-79</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.6e-77</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AL133412.1</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">5.2</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.40</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.0e-04</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.8e-03</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">HMGN2P17</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.7</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.48</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.2e-40</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">9.2e-39</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AC107959.5</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.6</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.25</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.8e-85</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.4e-82</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">CTBP2P6</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.6</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.65</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.1e-19</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.0e-18</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AC096861.2</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.5</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.50</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">5.8e-34</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.9e-32</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">PDLIM3</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.5</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.10</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.8e-04</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.3e-03</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">PPIAP74</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.1</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.14</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">9.0e-208</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">5.2e-205</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">CR383656.13</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.89</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.1e-05</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.0e-04</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">PPIAP8</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.8</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.18</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">5.8e-105</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.2e-102</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AC008391.1</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.7</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.70</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.5e-08</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.2e-06</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AC090735.1</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.4</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.33</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.3e-30</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.0e-29</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AC109129.1</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.4</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.49</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.9e-16</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.8e-14</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AL450998.1</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.4</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.52</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.7e-12</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.5e-11</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">CCDC192</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.3</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.50</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.4e-13</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">6.2e-12</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">MT1M</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.2</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.96</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">7.2e-03</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.1e-02</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">CFL1P2</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.1</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.35</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">7.1e-24</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.3e-22</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AMD1P4</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.35</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.4e-23</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">6.2e-22</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">PPIAP59</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.40</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.7e-17</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">5.9e-16</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">RN7SKP250</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.8</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.45</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">3.5e-12</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.1e-11</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6bba"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">AL139130.1</span></p></td><td class="cl-2c3f6ba6"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.7</span></p></td><td class="cl-2c3f6b6a"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.69</span></p></td><td class="cl-2c3f6bd8"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">8.0e-05</span></p></td><td class="cl-2c3f6bc4"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">6.3e-04</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-2c3f6c1e"><p class="cl-2c3ef7f2"><span class="cl-2c3edbe6">RPS4XP22</span></p></td><td class="cl-2c3f6bec"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">2.7</span></p></td><td class="cl-2c3f6c32"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">0.34</span></p></td><td class="cl-2c3f6c00"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">4.3e-17</span></p></td><td class="cl-2c3f6c14"><p class="cl-2c3ef824"><span class="cl-2c3edbe6">1.4e-15</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="5"class="cl-2c3f6c82"><p class="cl-2c3ef7f2"><span class="cl-2c3edc54">*</span><span class="cl-2c3edbe6">Wald test p-values</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="5"class="cl-2c3f6c82"><p class="cl-2c3ef7f2"><span class="cl-2c3edc54">†</span><span class="cl-2c3edbe6">Benjamini–Hochberg adjusted value</span></p></td></tr></tfoot></table></div></template>
<div id="e09feff3-4e27-4651-8b84-86f911ac1617"></div>
<script>
var dest = document.getElementById("e09feff3-4e27-4651-8b84-86f911ac1617");
var template = document.getElementById("0a79c07c-ba2b-4576-b32a-7f98b7256987");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;"
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
fantome.appendChild(templateContent);
</script>

```

\newpage

## Expression of top class DE genes by class {.tabset .tabset-fade .tabset-pills}



### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/top_class_DE_disease_sort-1.png" width="1152" />

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/top_class_DE_cluster_sort-1.png" width="1152" />

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/top_class_DE_hierarchical_sort-1.png" width="1152" />

\newpage

## Expression of all top DE genes {.tabset .tabset-fade .tabset-pills}




### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Expression_of_all_top_DE_genes_disease_sort-1.png" width="1152" />

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Expression_of_all_top_DE_genes_cluster_sort-1.png" width="1152" />

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Expression_of_all_top_DE_genes_auto_sort-1.png" width="1152" />

\newpage

## Volcano plot of gene expression
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Volcano_plot_of_gene_expression-1.png" width="864" />

\newpage

## Expression of selected interferon-stimulated genes {.tabset .tabset-fade .tabset-pills}


### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Expression_of_ISGs_disease_sort-1.png" width="1152" />


### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Expression_of_ISGs_cluster_sort-1.png" width="1152" />

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Expression_of_ISGs_auto_sort-1.png" width="1152" />

\newpage

## Expression of individual genes in the M1.2, M3.4, and M5.12 gene modules {.tabset .tabset-fade .tabset-pills}


### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/module_ISGs_disease_sort-1.png" width="1152" />

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/module_ISGs_cluster_sort-1.png" width="1152" />

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/module_ISGs_auto_sort-1.png" width="1152" />

\newpage

# Sample Clustering

The eigenvales for the annotated Banchereau modules were used to calculate an
adjacentcy matrix, which was in turn used to calculate the gap statistic to
determine the optimal *k*-clusters.

<!--html_preserve--><div id="htmlwidget-dfa74dc31f3af5d2f5ff" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-dfa74dc31f3af5d2f5ff">{"x":{"filter":"none","data":[["1","2"],["control","infected"],[16,12],[10,11],[23,14],[11,28],[16,14],[11,15]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>disease_class<\/th>\n      <th>1<\/th>\n      <th>2<\/th>\n      <th>3<\/th>\n      <th>4<\/th>\n      <th>5<\/th>\n      <th>6<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script><!--/html_preserve--><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Sample_Clustering_chart-2.png" width="672" />

\newpage

# Clustering and variation of samples in reduced dimensional space {.tabset .tabset-fade .tabset-pills}

## PCA {.tabset .tabset-fade .tabset-pills}


\newpage

### PCA by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/PCA_by_disease_classification-1.png" width="1152" />

\newpage

### PCA by plate
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/PCA_by_plate-1.png" width="1152" />

\newpage

### PCA by ethnicity
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/PCA_by_ethnicity-1.png" width="1152" />

\newpage

### PCA by K-means cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/PCA_by_K-means cluster-1.png" width="672" />

<!-- ### 3D PCA colored by K-means cluster -->
<!-- ```{r 3D_PCA_colored_by_K-means cluster} -->
<!-- pca_results %>% -->
<!--   plot_ly( -->
<!--     x = ~PC1, -->
<!--     y = ~PC2, -->
<!--     z = ~PC3, -->
<!--     colors = group_pal$cluster -->
<!--   ) %>% -->
<!--   add_markers( -->
<!--     text = str_glue("Sample: {pca_results$sample_name}<br>Classification: {pca_results$disease_class}<br>Ethnicity {pca_results$ethnicity}<br>Cluster: {pca_results$cluster}<br>Sex: {pca_results$sex}<br>Plate: {pca_results$run_id}"), -->
<!--     color = ~cluster, -->
<!--     size = 2.5, -->
<!--     alpha = 0.9, -->
<!--     hovertemplate = "%{text}", -->
<!--     marker = list( -->
<!--       size = 3, -->
<!--       line = list( -->
<!--         color = 'rgba(0, 0, 0, .8)', -->
<!--         width = 1.5 -->
<!--       ) -->
<!--     ) -->
<!--   ) -->
<!-- ``` -->

\newpage

## UMAP {.tabset .tabset-fade .tabset-pills}



### UMAP by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/UMAP_by_disease_classification-1.png" width="1152" />

\newpage

### UMAP by plate
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/UMAP_by_plate-1.png" width="1152" />

\newpage

### UMAP by ethnicity
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/UMAP_by_ethnicity-1.png" width="1152" />

\newpage

### UMAP by K-means cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/UMAP_by_K-means cluster-1.png" width="672" />

<!-- ### UMAP in 3D by K-means cluster -->
<!-- ```{r UMAP_in_3D_by_K-means cluster} -->
<!-- umap_results %>% -->
<!--   plot_ly( -->
<!--     x = ~umap_1, -->
<!--     y = ~umap_2, -->
<!--     z = ~umap_3, -->
<!--     colors = group_pal$cluster -->
<!--   ) %>% -->
<!--   add_markers( -->
<!--     text = str_glue("Sample: {umap_results$sample_name}<br>Classification: {umap_results$disease_class}<br>Ethnicity {umap_results$ethnicity}<br>Cluster: {umap_results$cluster}<br>Sex: {umap_results$sex}<br>Plate: {umap_results$run_id}"), -->
<!--     color = ~cluster, -->
<!--     size = 2.5, -->
<!--     alpha = 0.9, -->
<!--     hovertemplate = "%{text}", -->
<!--     marker = list( -->
<!--       size = 3, -->
<!--       line = list( -->
<!--         color = 'rgba(0, 0, 0, .8)', -->
<!--         width = 1.5 -->
<!--       ) -->
<!--     ) -->
<!--   ) -->
<!-- ``` -->

\newpage

# Module scores {.tabset .tabset-fade .tabset-pills}

## SLE Interferon module scores {.tabset .tabset-fade .tabset-pills}


### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/SLE_Interferon_module_scores_disease_sort-1.png" width="288" />

\newpage

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/SLE_Interferon_module_scores_cluster_sort-1.png" width="288" />

\newpage

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/SLE_Interferon_module_scores_auto_sort-1.png" width="384" />

\newpage

## SLE Inflammatory module scores {.tabset .tabset-fade .tabset-pills}


### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/SLE_Inflammatory_module_scores_disease_sort-1.png" width="384" />

\newpage

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/SLE_Inflammatory_module_scores_cluster_sort-1.png" width="384" />

\newpage

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/SLE_Inflammatory_module_scores_auto_sort-1.png" width="432" />

## Low-density granulocyte module scores {.tabset .tabset-fade .tabset-pills}


### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/ldg_module_disease_sort-1.png" width="312" />

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/ldg_module_cluster_sort-1.png" width="312" />

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/ldg_module_auto_sort-1.png" width="384" />

## SLE Modules with an annotated associated pathway {.tabset .tabset-fade .tabset-pills}





### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/annotated_module_disease_sort-1.png" width="1152" />

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/annotated_modulecluster_sort-1.png" width="1152" />

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/annotated_module_auto_sort-1.png" width="1152" />

\newpage

## All SLE modules {.tabset .tabset-fade .tabset-pills}


### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/all_module_disease_sort-1.png" width="1152" />

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/all_module_cluster_sort-1.png" width="1152" />

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/all_module_auto_sort-1.png" width="1152" />

\newpage

## Expression of the annotated Banchereau modules {.tabset .tabset-fade .tabset-pills}


### By disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_disease-1.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_disease-2.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_disease-3.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_disease-4.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_disease-5.png" width="1152" />

### By cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_cluster-1.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_cluster-2.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_cluster-3.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_cluster-4.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/Banchereau_modules_by_cluster-5.png" width="1152" />

\newpage

# WGCNA


## Expression similarity dendrogram
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/WGCNA_Expression_similarity_dendrogram-1.png" width="960" />

## Expression of identified modules {.tabset .tabset-fade .tabset-pills}


### By disease classification 
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_by_disease-1.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_by_disease-2.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_by_disease-3.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_by_disease-4.png" width="1152" />

\newpage

### By cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_modules_by_cluster-1.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_modules_by_cluster-2.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_modules_by_cluster-3.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_modules_by_cluster-4.png" width="1152" />

\newpage

## WGCNA eigenvalues {.tabset .tabset-fade .tabset-pills}


### Sorted by disease classification
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/WGCNA_eigenvalues_disease_sort-1.png" width="1152" />

### Sorted by cluster
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/WGCNA_eigenvalues_cluster_sort-1.png" width="1152" />

### Sorted by hierarchical clustering
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/WGCNA_eigenvalues_auto_sort-1.png" width="1152" />

\newpage

## Gene ontology GSEA of WGCNA modules
<img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_gsea-1.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_gsea-2.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_gsea-3.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_gsea-4.png" width="1152" /><img src="/home/rstudio/oldworkspace/analysis/rnaseq/covid/covid_pcv/reports/osctr_report_files/figure-html/wgcna_gsea-5.png" width="1152" />

\newpage

## Use of WGCNA modules in classifying clusters.

The module eigengene scores were used to train a random forest model with repeated cross validation to predict cluster identity.

The model was trained with 75% of the dataset and tested on 25%.

Note: Subjects with LP were excluded due to the small population size.


```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction 1 2 3 4 5 6
##          1 3 0 1 1 0 0
##          2 0 2 0 0 0 0
##          3 0 3 6 0 0 3
##          4 3 0 0 3 0 0
##          5 0 0 0 5 7 0
##          6 1 0 2 0 0 3
## 
## Overall Statistics
##                                           
##                Accuracy : 0.5581          
##                  95% CI : (0.3988, 0.7092)
##     No Information Rate : 0.2093          
##     P-Value [Acc > NIR] : 5.746e-07       
##                                           
##                   Kappa : 0.4632          
##                                           
##  Mcnemar's Test P-Value : NA              
## 
## Statistics by Class:
## 
##                      Class: 1 Class: 2 Class: 3 Class: 4 Class: 5 Class: 6
## Sensitivity           0.42857  0.40000   0.6667  0.33333   1.0000  0.50000
## Specificity           0.94444  1.00000   0.8235  0.91176   0.8611  0.91892
## Pos Pred Value        0.60000  1.00000   0.5000  0.50000   0.5833  0.50000
## Neg Pred Value        0.89474  0.92683   0.9032  0.83784   1.0000  0.91892
## Prevalence            0.16279  0.11628   0.2093  0.20930   0.1628  0.13953
## Detection Rate        0.06977  0.04651   0.1395  0.06977   0.1628  0.06977
## Detection Prevalence  0.11628  0.04651   0.2791  0.13953   0.2791  0.13953
## Balanced Accuracy     0.68651  0.70000   0.7451  0.62255   0.9306  0.70946
```

The relative importance of each eigengene in classification:


```
## parRF variable importance
## 
##   variables are sorted by maximum importance across the classes
##   only 20 most important variables shown (out of 29)
## 
##                        1       2        3       4       5        6
## MEgreen          5.57200  9.9289  8.83706  4.3115 12.7264  3.44071
## MEpurple         0.73291  7.7817  4.37087  1.3100 12.6255  1.94943
## MEyellow         3.20920 10.5486  8.74307  9.2054 12.3572  4.73450
## MEcyan           2.60041  8.1296  1.89282  5.7920  9.6695  0.30893
## MEgreenyellow    0.59255  5.9454  2.89713  8.6253  4.9131  1.71677
## MEblack          0.21738  0.2292  2.11458  8.4522 -1.4980  0.72657
## MEmagenta        0.80775  4.1830  4.60225  7.3760  1.5807  1.29998
## MEdarkgreen      0.67047  5.6226  4.30671  2.4206  6.7844  1.57468
## MEwhite          0.45935  5.4996  0.05911  3.6010  3.1108 -0.25605
## MEgrey60        -0.48700  5.2622  1.33757  1.6807  1.5968  2.37445
## MEblue           0.88190  5.1609  1.97210  4.1214  4.9512  0.63125
## MEpink           3.07331  2.0035  1.56049  4.9580 -0.2291  2.67851
## MEturquoise     -1.09825 -0.1622  2.69052  4.3309 -0.9887 -1.97860
## MEdarkred       -1.00735  3.5835 -0.12616  0.2841  4.3264 -2.02659
## MEsalmon        -0.06547  4.2067 -0.09559  0.4372  1.0076 -0.57672
## MElightgreen    -1.89724  1.3325  3.21842 -0.8997  1.9923  4.12991
## MEtan            0.42020  1.1131 -0.44010  3.7725  4.1108  0.99891
## MEdarkturquoise  0.53502 -0.2410 -0.15567  2.7985  3.8976  0.58223
## MElightyellow    1.82939  3.1575  0.96367  3.2648 -0.3915  1.12035
## MEmidnightblue   1.80007  2.9268  2.94486  2.6207  3.0821 -0.07906
```

\newpage

## Use of WGCNA modules in labeling disease class.

The module eigengene scores were used to train a random forest model with repeated cross validation to predict cluster identity.

The model was trained with 75% of the dataset and tested on 25%.


```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction control infected
##   control       21        0
##   infected       0       23
##                                      
##                Accuracy : 1          
##                  95% CI : (0.9196, 1)
##     No Information Rate : 0.5227     
##     P-Value [Acc > NIR] : 4.019e-13  
##                                      
##                   Kappa : 1          
##                                      
##  Mcnemar's Test P-Value : NA         
##                                      
##             Sensitivity : 1.0000     
##             Specificity : 1.0000     
##          Pos Pred Value : 1.0000     
##          Neg Pred Value : 1.0000     
##              Prevalence : 0.4773     
##          Detection Rate : 0.4773     
##    Detection Prevalence : 0.4773     
##       Balanced Accuracy : 1.0000     
##                                      
##        'Positive' Class : control    
## 
```

The relative importance of each eigengene in classification:


```
## parRF variable importance
## 
##   only 20 most important variables shown (out of 29)
## 
##                Overall
## MEturquoise    48.1414
## MEbrown         6.9179
## MEmagenta       5.0464
## MEblue          1.4670
## MEpink          0.9205
## MEblack         0.5741
## MEmidnightblue  0.5027
## MEsalmon        0.4836
## MEskyblue       0.3459
## MEdarkred       0.3280
## MEsaddlebrown   0.3149
## MEwhite         0.3011
## MEroyalblue     0.2742
## MEcyan          0.2602
## MEyellow        0.2461
## MElightgreen    0.2216
## MEtan           0.1763
## MEdarkgrey      0.1669
## MEpurple        0.1645
## MEgreenyellow   0.1480
```

\newpage

## Use of module modules in classifying clusters.

The module eigengene scores were used to train a random forest model with repeated cross validation to predict cluster identity.

The model was trained with 75% of the dataset and tested on 25%.


```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction 1 2 3 4 5 6
##          1 4 0 0 0 0 1
##          2 0 4 0 0 0 0
##          3 0 1 9 0 0 1
##          4 1 0 0 8 2 0
##          5 0 0 0 1 5 0
##          6 2 0 0 0 0 4
## 
## Overall Statistics
##                                           
##                Accuracy : 0.7907          
##                  95% CI : (0.6396, 0.8996)
##     No Information Rate : 0.2093          
##     P-Value [Acc > NIR] : 5.888e-16       
##                                           
##                   Kappa : 0.7451          
##                                           
##  Mcnemar's Test P-Value : NA              
## 
## Statistics by Class:
## 
##                      Class: 1 Class: 2 Class: 3 Class: 4 Class: 5 Class: 6
## Sensitivity           0.57143  0.80000   1.0000   0.8889   0.7143  0.66667
## Specificity           0.97222  1.00000   0.9412   0.9118   0.9722  0.94595
## Pos Pred Value        0.80000  1.00000   0.8182   0.7273   0.8333  0.66667
## Neg Pred Value        0.92105  0.97436   1.0000   0.9687   0.9459  0.94595
## Prevalence            0.16279  0.11628   0.2093   0.2093   0.1628  0.13953
## Detection Rate        0.09302  0.09302   0.2093   0.1860   0.1163  0.09302
## Detection Prevalence  0.11628  0.09302   0.2558   0.2558   0.1395  0.13953
## Balanced Accuracy     0.77183  0.90000   0.9706   0.9003   0.8433  0.80631
```

The relative importance of each eigengene in classification:


```
## parRF variable importance
## 
##   variables are sorted by maximum importance across the classes
##   only 20 most important variables shown (out of 260)
## 
##              1     2        3       4     5        6
## M6.13  3.18171 2.996  3.50258  2.3137 4.137  0.92408
## M8.2   0.18517 4.128 -0.03508  1.5046 3.201  2.69942
## M7.1   1.69880 4.095  1.28358  3.4773 3.348  3.34142
## M7.15  0.07025 3.974  2.52504  2.7136 3.997  1.14740
## M5.7   0.34483 2.827  1.73688  3.1309 3.984  0.62653
## M5.14  0.71304 3.888  3.94841  2.2197 3.598 -1.27518
## M3.2   1.40115 3.887  3.03085  1.9328 3.825  0.88991
## M8.10  0.76884 3.742  2.49105  2.6819 3.843  1.76012
## M8.48 -0.62201 3.052  1.08804  1.9219 3.813  0.72382
## M7.12 -0.24227 3.799  1.00117  2.4786 2.606  2.17807
## M7.27  1.23095 3.764  3.78302  2.8899 2.685  1.07033
## M4.9   1.04093 1.104  3.26132  1.9804 3.765  0.08432
## M8.88  1.93657 3.705  2.09763  2.1965 3.129 -0.38453
## M4.13  1.65464 3.483  3.28441  2.5049 3.672  0.04231
## M6.12  0.09992 3.444  3.24888  1.2090 3.669  0.91065
## M8.14  2.49561 3.293  0.11507 -0.4465 3.640  1.14549
## M9.15  1.50053 2.548  2.53657  0.4105 3.621  0.96827
## M9.22  0.43153 2.687  2.03901 -1.4961 3.573  1.46419
## M4.2   0.85894 2.944  1.17188  0.9191 3.557  1.84895
## M4.6   2.10532 3.546  2.72429  2.6930 3.529  1.95350
```

\newpage

## Use of module modules in labeling disease class.

The module eigengene scores were used to train a random forest model with repeated cross validation to predict cluster identity.

The model was trained with 75% of the dataset and tested on 25%.


```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction control infected
##   control       17        2
##   infected       4       21
##                                           
##                Accuracy : 0.8636          
##                  95% CI : (0.7265, 0.9483)
##     No Information Rate : 0.5227          
##     P-Value [Acc > NIR] : 1.963e-06       
##                                           
##                   Kappa : 0.7256          
##                                           
##  Mcnemar's Test P-Value : 0.6831          
##                                           
##             Sensitivity : 0.8095          
##             Specificity : 0.9130          
##          Pos Pred Value : 0.8947          
##          Neg Pred Value : 0.8400          
##              Prevalence : 0.4773          
##          Detection Rate : 0.3864          
##    Detection Prevalence : 0.4318          
##       Balanced Accuracy : 0.8613          
##                                           
##        'Positive' Class : control         
## 
```

The relative importance of each eigengene in classification:


```
## parRF variable importance
## 
##   only 20 most important variables shown (out of 260)
## 
##        Importance
## M8.50      16.428
## M8.21      11.281
## M6.2        7.138
## M7.19       6.076
## M7.7        4.758
## M7.8        4.473
## M8.36       4.148
## M9.7        3.711
## M5.14       3.685
## M7.18       3.274
## M4.11       3.026
## M5.6        2.786
## M6.5        2.711
## M8.95       2.707
## M5.10       2.687
## M6.16       2.616
## M8.107      2.434
## M9.37       2.380
## M9.26       2.304
## M8.66       2.285
```


```
## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 4.0.3 (2020-10-10)
##  os       Ubuntu 20.04.1 LTS          
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  ctype    en_US.UTF-8                 
##  tz       Etc/UTC                     
##  date     2021-03-15                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  package              * version    date       lib source                                   
##  abind                  1.4-5      2016-07-21 [1] RSPM (R 4.0.0)                           
##  annotate               1.68.0     2020-10-27 [1] Bioconductor                             
##  AnnotationDbi          1.52.0     2020-10-27 [1] Bioconductor                             
##  assertthat             0.2.1      2019-03-21 [1] RSPM (R 4.0.0)                           
##  backports              1.2.1      2020-12-09 [1] RSPM (R 4.0.3)                           
##  base64enc              0.1-3      2015-07-28 [1] RSPM (R 4.0.0)                           
##  base64url              1.4        2018-05-14 [1] RSPM (R 4.0.0)                           
##  Biobase              * 2.50.0     2020-10-27 [1] Bioconductor                             
##  BiocGenerics         * 0.36.0     2020-10-27 [1] Bioconductor                             
##  BiocManager          * 1.30.10    2019-11-16 [1] CRAN (R 4.0.3)                           
##  BiocParallel         * 1.24.1     2020-11-06 [1] Bioconductor                             
##  bit                    4.0.4      2020-08-04 [1] RSPM (R 4.0.2)                           
##  bit64                  4.0.5      2020-08-30 [1] RSPM (R 4.0.2)                           
##  bitops                 1.0-6      2013-08-17 [1] RSPM (R 4.0.0)                           
##  blob                   1.2.1      2020-01-20 [1] RSPM (R 4.0.0)                           
##  broom                * 0.7.3      2020-12-16 [1] RSPM (R 4.0.3)                           
##  car                    3.0-10     2020-09-29 [1] RSPM (R 4.0.2)                           
##  carData                3.0-4      2020-05-22 [1] RSPM (R 4.0.0)                           
##  caret                * 6.0-86     2020-03-20 [1] RSPM (R 4.0.0)                           
##  cellranger             1.1.0      2016-07-27 [1] RSPM (R 4.0.0)                           
##  checkmate              2.0.0      2020-02-06 [1] RSPM (R 4.0.0)                           
##  class                  7.3-17     2020-04-26 [2] CRAN (R 4.0.3)                           
##  cli                    2.3.1      2021-02-23 [1] RSPM (R 4.0.3)                           
##  cluster              * 2.1.0      2019-06-19 [1] RSPM (R 4.0.0)                           
##  clusterProfiler      * 3.18.0     2020-10-27 [1] Bioconductor                             
##  codetools              0.2-18     2020-11-04 [2] RSPM (R 4.0.3)                           
##  colorspace             2.0-0      2020-11-11 [1] RSPM (R 4.0.3)                           
##  corrplot             * 0.84       2017-10-16 [1] RSPM (R 4.0.0)                           
##  cowplot              * 1.1.1      2020-12-30 [1] RSPM (R 4.0.3)                           
##  crayon                 1.4.1      2021-02-08 [1] RSPM (R 4.0.3)                           
##  crosstalk              1.1.0.1    2020-03-13 [1] RSPM (R 4.0.0)                           
##  curl                   4.3        2019-12-02 [1] RSPM (R 4.0.0)                           
##  data.table           * 1.13.6     2020-12-30 [1] RSPM (R 4.0.3)                           
##  DBI                    1.1.0      2019-12-15 [1] RSPM (R 4.0.0)                           
##  dbplyr                 2.0.0      2020-11-03 [1] RSPM (R 4.0.3)                           
##  DelayedArray           0.16.0     2020-10-27 [1] Bioconductor                             
##  DESeq2               * 1.30.0     2020-10-27 [1] Bioconductor                             
##  dials                * 0.0.9      2020-09-16 [1] RSPM (R 4.0.2)                           
##  DiceDesign             1.8-1      2019-07-31 [1] RSPM (R 4.0.0)                           
##  digest                 0.6.27     2020-10-24 [1] RSPM (R 4.0.3)                           
##  DO.db                  2.9        2021-01-02 [1] Bioconductor                             
##  doParallel             1.0.16     2020-10-16 [1] RSPM (R 4.0.2)                           
##  DOSE                   3.16.0     2020-10-27 [1] Bioconductor                             
##  downloader             0.4        2015-07-09 [1] RSPM (R 4.0.0)                           
##  dplyr                * 1.0.2      2020-08-18 [1] RSPM (R 4.0.2)                           
##  drake                * 7.12.7     2020-10-27 [1] RSPM (R 4.0.3)                           
##  DT                     0.16       2020-10-13 [1] RSPM (R 4.0.2)                           
##  dynamicTreeCut       * 1.63-1     2016-03-11 [1] RSPM (R 4.0.0)                           
##  e1071                  1.7-4      2020-10-14 [1] RSPM (R 4.0.3)                           
##  edgeR                  3.32.0     2020-10-27 [1] Bioconductor                             
##  ellipsis               0.3.1      2020-05-15 [1] RSPM (R 4.0.0)                           
##  enrichplot             1.10.1     2020-11-14 [1] Bioconductor                             
##  evaluate               0.14       2019-05-28 [1] RSPM (R 4.0.0)                           
##  factoextra           * 1.0.7      2020-04-01 [1] RSPM (R 4.0.0)                           
##  fansi                  0.4.2      2021-01-15 [1] RSPM (R 4.0.3)                           
##  farver                 2.1.0      2021-02-28 [1] RSPM (R 4.0.3)                           
##  fastcluster          * 1.1.25     2018-06-07 [1] RSPM (R 4.0.0)                           
##  fastmatch              1.1-0      2017-01-28 [1] RSPM (R 4.0.0)                           
##  fgsea                  1.16.0     2020-10-27 [1] Bioconductor                             
##  filelock               1.0.2      2018-10-05 [1] RSPM (R 4.0.0)                           
##  flextable            * 0.6.1      2020-12-08 [1] RSPM (R 4.0.3)                           
##  forcats              * 0.5.0      2020-03-01 [1] RSPM (R 4.0.0)                           
##  foreach                1.5.1      2020-10-15 [1] RSPM (R 4.0.2)                           
##  foreign                0.8-81     2020-12-22 [2] RSPM (R 4.0.3)                           
##  formattable          * 0.2.0.1    2016-08-05 [1] RSPM (R 4.0.0)                           
##  Formula                1.2-4      2020-10-16 [1] RSPM (R 4.0.2)                           
##  fs                     1.5.0      2020-07-31 [1] RSPM (R 4.0.2)                           
##  furrr                * 0.2.1      2020-10-21 [1] RSPM (R 4.0.3)                           
##  future               * 1.21.0     2020-12-10 [1] RSPM (R 4.0.3)                           
##  gdtools                0.2.2      2020-04-03 [1] RSPM (R 4.0.2)                           
##  genefilter           * 1.72.0     2020-10-27 [1] Bioconductor                             
##  geneplotter            1.68.0     2020-10-27 [1] Bioconductor                             
##  generics               0.1.0      2020-10-31 [1] RSPM (R 4.0.3)                           
##  GenomeInfoDb         * 1.26.2     2020-12-08 [1] Bioconductor                             
##  GenomeInfoDbData       1.2.4      2021-01-02 [1] Bioconductor                             
##  GenomicRanges        * 1.42.0     2020-10-27 [1] Bioconductor                             
##  ggfittext              0.9.0      2020-06-14 [1] RSPM (R 4.0.0)                           
##  ggforce              * 0.3.2      2020-06-23 [1] RSPM (R 4.0.3)                           
##  ggplot2              * 3.3.3      2020-12-30 [1] RSPM (R 4.0.3)                           
##  ggplotify            * 0.0.5      2020-03-12 [1] RSPM (R 4.0.0)                           
##  ggpubr               * 0.4.0      2020-06-27 [1] RSPM (R 4.0.2)                           
##  ggradar              * 0.2        2021-01-02 [1] Github (ricardo-bion/ggradar@63e5cef)    
##  ggraph                 2.0.4      2020-11-16 [1] RSPM (R 4.0.3)                           
##  ggrepel              * 0.9.0      2020-12-16 [1] RSPM (R 4.0.3)                           
##  ggsignif               0.6.0      2019-08-08 [1] RSPM (R 4.0.0)                           
##  ggtext               * 0.1.1      2020-12-17 [1] RSPM (R 4.0.3)                           
##  globals                0.14.0     2020-11-22 [1] RSPM (R 4.0.3)                           
##  glue                   1.4.2      2020-08-27 [1] RSPM (R 4.0.2)                           
##  GO.db                  3.12.1     2021-01-02 [1] Bioconductor                             
##  GOSemSim               2.16.1     2020-10-29 [1] Bioconductor                             
##  gower                  0.2.2      2020-06-23 [1] RSPM (R 4.0.2)                           
##  GPfit                  1.0-8      2019-02-08 [1] RSPM (R 4.0.0)                           
##  graphlayouts           0.7.1      2020-10-26 [1] RSPM (R 4.0.3)                           
##  gridBase               0.4-7      2014-02-24 [1] RSPM (R 4.0.0)                           
##  gridExtra              2.3        2017-09-09 [1] RSPM (R 4.0.0)                           
##  gridGraphics           0.5-1      2020-12-13 [1] RSPM (R 4.0.3)                           
##  gridtext               0.1.4      2020-12-10 [1] RSPM (R 4.0.3)                           
##  gtable                 0.3.0      2019-03-25 [1] RSPM (R 4.0.0)                           
##  gtools               * 3.8.2      2020-03-31 [1] RSPM (R 4.0.0)                           
##  haven                  2.3.1      2020-06-01 [1] RSPM (R 4.0.2)                           
##  here                 * 1.0.1      2020-12-13 [1] RSPM (R 4.0.3)                           
##  HGNChelper           * 0.8.1      2019-10-24 [1] RSPM (R 4.0.0)                           
##  Hmisc                  4.4-2      2020-11-29 [1] RSPM (R 4.0.3)                           
##  hms                    0.5.3      2020-01-08 [1] RSPM (R 4.0.0)                           
##  htmlTable              2.1.0      2020-09-16 [1] RSPM (R 4.0.2)                           
##  htmltools              0.5.0      2020-06-16 [1] RSPM (R 4.0.1)                           
##  htmlwidgets            1.5.3      2020-12-10 [1] RSPM (R 4.0.3)                           
##  httr                   1.4.2      2020-07-20 [1] RSPM (R 4.0.2)                           
##  igraph                 1.2.6      2020-10-06 [1] RSPM (R 4.0.2)                           
##  import                 1.2.0      2020-09-24 [1] RSPM (R 4.0.2)                           
##  impute                 1.64.0     2020-10-27 [1] Bioconductor                             
##  infer                * 0.5.3      2020-07-14 [1] RSPM (R 4.0.2)                           
##  inspectdf            * 0.0.9      2020-09-07 [1] RSPM (R 4.0.2)                           
##  ipred                  0.9-9      2019-04-28 [1] RSPM (R 4.0.0)                           
##  IRanges              * 2.24.1     2020-12-12 [1] Bioconductor                             
##  irlba                * 2.3.3      2019-02-05 [1] RSPM (R 4.0.3)                           
##  iterators              1.0.13     2020-10-15 [1] RSPM (R 4.0.2)                           
##  janitor              * 2.0.1      2020-04-12 [1] RSPM (R 4.0.0)                           
##  jpeg                   0.1-8.1    2019-10-24 [1] RSPM (R 4.0.0)                           
##  jsonlite               1.7.2      2020-12-09 [1] RSPM (R 4.0.3)                           
##  kableExtra           * 1.3.1      2020-10-22 [1] RSPM (R 4.0.3)                           
##  knitr                * 1.30       2020-09-22 [1] RSPM (R 4.0.2)                           
##  labeling               0.4.2      2020-10-20 [1] RSPM (R 4.0.3)                           
##  lattice              * 0.20-41    2020-04-02 [2] CRAN (R 4.0.3)                           
##  latticeExtra           0.6-29     2019-12-19 [1] RSPM (R 4.0.0)                           
##  lava                   1.6.8.1    2020-11-04 [1] RSPM (R 4.0.3)                           
##  lazyeval               0.2.2      2019-03-15 [1] RSPM (R 4.0.0)                           
##  lhs                    1.1.1      2020-10-05 [1] RSPM (R 4.0.2)                           
##  lifecycle              1.0.0      2021-02-15 [1] RSPM (R 4.0.3)                           
##  limma                  3.46.0     2020-10-27 [1] Bioconductor                             
##  listenv                0.8.0      2019-12-05 [1] RSPM (R 4.0.0)                           
##  locfit                 1.5-9.4    2020-03-25 [1] RSPM (R 4.0.0)                           
##  lubridate              1.7.9.2    2020-11-13 [1] RSPM (R 4.0.3)                           
##  magrittr             * 2.0.1      2020-11-17 [1] RSPM (R 4.0.3)                           
##  markdown               1.1        2019-08-07 [1] RSPM (R 4.0.0)                           
##  MASS                   7.3-53     2020-09-09 [2] CRAN (R 4.0.3)                           
##  Matrix               * 1.3-0      2020-12-22 [2] RSPM (R 4.0.3)                           
##  MatrixGenerics       * 1.2.0      2020-10-27 [1] Bioconductor                             
##  matrixStats          * 0.57.0     2020-09-25 [1] RSPM (R 4.0.2)                           
##  memoise                1.1.0      2017-04-21 [1] RSPM (R 4.0.0)                           
##  mgcv                 * 1.8-33     2020-08-27 [2] CRAN (R 4.0.3)                           
##  modeldata            * 0.1.0      2020-10-22 [1] RSPM (R 4.0.3)                           
##  ModelMetrics           1.2.2.2    2020-03-17 [1] RSPM (R 4.0.2)                           
##  modelr                 0.1.8      2020-05-19 [1] RSPM (R 4.0.0)                           
##  moduleScoreR         * 0.0.0.9200 2021-01-02 [1] Github (milescsmith/moduleScoreR@e6db70f)
##  munsell                0.5.0      2018-06-12 [1] RSPM (R 4.0.0)                           
##  nlme                 * 3.1-151    2020-12-10 [2] RSPM (R 4.0.3)                           
##  NMF                    0.23.0     2020-08-01 [1] RSPM (R 4.0.3)                           
##  nnet                   7.3-14     2020-04-26 [2] CRAN (R 4.0.3)                           
##  oaColors             * 0.0.4      2015-11-30 [1] RSPM (R 4.0.0)                           
##  officer                0.3.15     2020-11-01 [1] RSPM (R 4.0.3)                           
##  openxlsx               4.2.3      2020-10-27 [1] RSPM (R 4.0.3)                           
##  paletteer            * 1.2.0      2020-06-07 [1] RSPM (R 4.0.0)                           
##  parallelly             1.22.0     2020-12-13 [1] RSPM (R 4.0.3)                           
##  parsnip              * 0.1.4      2020-10-27 [1] RSPM (R 4.0.3)                           
##  pheatmap             * 1.0.12     2019-01-04 [1] RSPM (R 4.0.0)                           
##  pillar                 1.5.1      2021-03-05 [1] RSPM (R 4.0.3)                           
##  pkgconfig              2.0.3      2019-09-22 [1] RSPM (R 4.0.0)                           
##  pkgmaker               0.32.2     2020-10-20 [1] RSPM (R 4.0.3)                           
##  plotly               * 4.9.2.2    2020-12-19 [1] RSPM (R 4.0.3)                           
##  plyr                   1.8.6      2020-03-03 [1] RSPM (R 4.0.2)                           
##  png                    0.1-7      2013-12-03 [1] RSPM (R 4.0.0)                           
##  polyclip               1.10-0     2019-03-14 [1] RSPM (R 4.0.0)                           
##  preprocessCore         1.52.0     2020-10-27 [1] Bioconductor                             
##  prettyunits            1.1.1      2020-01-24 [1] RSPM (R 4.0.0)                           
##  prismatic              0.2.0      2019-12-01 [1] RSPM (R 4.0.0)                           
##  pROC                   1.16.2     2020-03-19 [1] RSPM (R 4.0.2)                           
##  prodlim                2019.11.13 2019-11-17 [1] RSPM (R 4.0.2)                           
##  progress               1.2.2      2019-05-16 [1] RSPM (R 4.0.0)                           
##  purrr                * 0.3.4      2020-04-17 [1] RSPM (R 4.0.0)                           
##  qvalue                 2.22.0     2020-10-27 [1] Bioconductor                             
##  R6                     2.5.0      2020-10-28 [1] RSPM (R 4.0.3)                           
##  randomForest         * 4.6-14     2018-03-25 [1] RSPM (R 4.0.0)                           
##  RColorBrewer         * 1.1-2      2014-12-07 [1] RSPM (R 4.0.0)                           
##  Rcpp                   1.0.6      2021-01-15 [1] RSPM (R 4.0.3)                           
##  RCurl                  1.98-1.2   2020-04-18 [1] RSPM (R 4.0.0)                           
##  readr                * 1.4.0      2020-10-05 [1] RSPM (R 4.0.3)                           
##  readxl               * 1.3.1      2019-03-13 [1] RSPM (R 4.0.2)                           
##  recipes              * 0.1.15     2020-11-11 [1] RSPM (R 4.0.3)                           
##  registry               0.5-1      2019-03-05 [1] RSPM (R 4.0.0)                           
##  rematch2               2.1.2      2020-05-01 [1] RSPM (R 4.0.0)                           
##  reprex                 0.3.0      2019-05-16 [1] RSPM (R 4.0.0)                           
##  reshape2               1.4.4      2020-04-09 [1] RSPM (R 4.0.2)                           
##  rio                    0.5.16     2018-11-26 [1] RSPM (R 4.0.0)                           
##  rlang                * 0.4.10     2020-12-30 [1] RSPM (R 4.0.3)                           
##  rmarkdown            * 2.6        2020-12-14 [1] RSPM (R 4.0.3)                           
##  rngtools               1.5        2020-01-23 [1] RSPM (R 4.0.0)                           
##  rpart                  4.1-15     2019-04-12 [2] CRAN (R 4.0.3)                           
##  rprojroot              2.0.2      2020-11-15 [1] RSPM (R 4.0.3)                           
##  rsample              * 0.0.8      2020-09-23 [1] RSPM (R 4.0.2)                           
##  RSQLite                2.2.1      2020-09-30 [1] RSPM (R 4.0.2)                           
##  rstatix              * 0.6.0      2020-06-18 [1] RSPM (R 4.0.1)                           
##  rstudioapi             0.13       2020-11-12 [1] RSPM (R 4.0.3)                           
##  rsvd                   1.0.3      2020-02-17 [1] RSPM (R 4.0.0)                           
##  rvcheck                0.1.8      2020-03-01 [1] RSPM (R 4.0.0)                           
##  rvest                  0.3.6      2020-07-25 [1] RSPM (R 4.0.2)                           
##  S4Vectors            * 0.28.1     2020-12-09 [1] Bioconductor                             
##  scales               * 1.1.1      2020-05-11 [1] RSPM (R 4.0.0)                           
##  scatterpie             0.1.5      2020-09-09 [1] RSPM (R 4.0.2)                           
##  sessioninfo            1.1.1      2018-11-05 [1] RSPM (R 4.0.0)                           
##  shadowtext             0.0.7      2019-11-06 [1] RSPM (R 4.0.0)                           
##  snakecase              0.11.0     2019-05-25 [1] RSPM (R 4.0.0)                           
##  snow                   0.4-3      2018-09-14 [1] RSPM (R 4.0.0)                           
##  storr                  1.2.5      2020-12-01 [1] RSPM (R 4.0.3)                           
##  stringi                1.5.3      2020-09-09 [1] RSPM (R 4.0.2)                           
##  stringr              * 1.4.0      2019-02-10 [1] RSPM (R 4.0.0)                           
##  SummarizedExperiment * 1.20.0     2020-10-27 [1] Bioconductor                             
##  survival               3.2-7      2020-09-28 [2] CRAN (R 4.0.3)                           
##  sva                  * 3.38.0     2020-10-27 [1] Bioconductor                             
##  systemfonts            0.3.2      2020-09-29 [1] RSPM (R 4.0.3)                           
##  tibble               * 3.1.0      2021-02-25 [1] RSPM (R 4.0.3)                           
##  tidygraph              1.2.0      2020-05-12 [1] RSPM (R 4.0.2)                           
##  tidymodels           * 0.1.2      2020-11-22 [1] RSPM (R 4.0.3)                           
##  tidyr                * 1.1.2      2020-08-27 [1] RSPM (R 4.0.3)                           
##  tidyselect             1.1.0      2020-05-11 [1] RSPM (R 4.0.0)                           
##  tidyverse            * 1.3.0      2019-11-21 [1] RSPM (R 4.0.0)                           
##  timeDate               3043.102   2018-02-21 [1] RSPM (R 4.0.0)                           
##  tune                 * 0.1.2      2020-11-17 [1] RSPM (R 4.0.3)                           
##  tweenr                 1.0.1      2018-12-14 [1] RSPM (R 4.0.2)                           
##  tximport             * 1.18.0     2020-10-27 [1] Bioconductor                             
##  txtq                   0.2.3      2020-06-23 [1] RSPM (R 4.0.2)                           
##  utf8                   1.1.4      2018-05-24 [1] RSPM (R 4.0.0)                           
##  uuid                   0.1-4      2020-02-26 [1] RSPM (R 4.0.0)                           
##  uwot                 * 0.1.10     2020-12-15 [1] RSPM (R 4.0.3)                           
##  vctrs                  0.3.6      2020-12-17 [1] RSPM (R 4.0.3)                           
##  viridis              * 0.5.1      2018-03-29 [1] RSPM (R 4.0.0)                           
##  viridisLite          * 0.3.0      2018-02-01 [1] RSPM (R 4.0.0)                           
##  webshot                0.5.2      2019-11-22 [1] RSPM (R 4.0.0)                           
##  WGCNA                * 1.69       2020-02-28 [1] RSPM (R 4.0.3)                           
##  withr                  2.4.1      2021-01-26 [1] RSPM (R 4.0.3)                           
##  workflows            * 0.2.1      2020-10-08 [1] RSPM (R 4.0.2)                           
##  xfun                   0.19       2020-10-30 [1] RSPM (R 4.0.3)                           
##  XML                    3.99-0.5   2020-07-23 [1] RSPM (R 4.0.2)                           
##  xml2                   1.3.2      2020-04-23 [1] RSPM (R 4.0.0)                           
##  xtable                 1.8-4      2019-04-21 [1] RSPM (R 4.0.0)                           
##  XVector                0.30.0     2020-10-27 [1] Bioconductor                             
##  yaml                   2.2.1      2020-02-01 [1] RSPM (R 4.0.0)                           
##  yardstick            * 0.0.7      2020-07-13 [1] RSPM (R 4.0.2)                           
##  zip                    2.1.1      2020-08-27 [1] RSPM (R 4.0.2)                           
##  zlibbioc               1.36.0     2020-10-27 [1] Bioconductor                             
## 
## [1] /usr/local/lib/R/site-library
## [2] /usr/local/lib/R/library
```
