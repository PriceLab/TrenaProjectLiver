library(TrenaViz)
tv <- TrenaViz("TrenaProjectLiver")
later(function(){browseURL("http://0.0.0.0:6838")}, 5)
runApp(createApp(tv, port=6838))
