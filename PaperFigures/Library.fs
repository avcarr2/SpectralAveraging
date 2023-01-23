namespace PaperFigures
open Tests
open Plotly.NET; 
open Plotly.NET.GenericChart
open Plotly.NET.LayoutObjects

module public PlotModifiers = 
    
    let ChangeAxesName (chart: GenericChart) xAxisTitle yAxisTitle : GenericChart = 
            let xAxis = 
                let tmp = LinearAxis()
                tmp?title <- xAxisTitle
                tmp?showGrid <- false
                tmp?showLine <- true
                tmp

            let yAxis = 
                let tmp = LinearAxis()
                tmp?title <- yAxisTitle
                tmp?showgrid <- false
                tmp?showline <- true
                tmp
            let layout = 
                let temp = Layout()
                temp?xaxis <- xAxis
                temp?yaxis <- yAxis
                temp?showLegend <- false
                temp
            chart
            |> GenericChart.addLayout layout 

        