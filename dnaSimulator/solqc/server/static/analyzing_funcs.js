
var interval;
$(document).ready(function () {

        $.blockUI({ message: $("#processing_div")});
        var time_stamp = document.getElementById("time_stamp").innerHTML ;
        interval = setInterval(()=> {
                var response =  $.get('/processing/' + time_stamp.toString(),
                    (data)=> {
                        console.log(data);
                        if(data === 'still working'){
                        //    do nothing maybe change to better solution
                        }
                        else if( data === 'error')
                        {
                            console.log(data);
                            clearInterval(interval);
                            var newDoc = document.open("text/html", "replace");
                            newDoc.write('Sorry, something went wrong');
                            newDoc.close();

                        }
                        //in case everything went properly (hopefully)
                        else {
                            console.log("i am so happy now")
                            console.log(data);
                            clearInterval(interval);
                            $.unblockUI();
                            var newDoc = document.open("text/html", "replace");
                            newDoc.write(data);
                            newDoc.close();
                        }

                    });


        }, 2000);

    });
