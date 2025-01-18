//this function checks all checkbox.
treat_as_radio_button = ['design_file_server', 'config_file_server']
mapp_row_name_to_value = {"design_file": {
                                            "url":"design_file_url",
                                            "upload":"design_file_upload",
                                            "mounted":"design_file_server"},
                            "config_file":{
                                            "url":"config_file_url",
                                            "upload":"config_file_upload",
                                            "mounted":"config_file_server"
                            },
                            "fastq_file":{
                                            "url":"zip_file_url",
                                            "upload":"zip_file_upload",
                                            "mounted":"fastq_file_server"
                            }

    };
function toggle(source) {
    console.log('Inside toggle');
  checkboxes = document.getElementsByName('fastq_file_server');

  console.log('checkboxes:',checkboxes)
  for(var i=0; i < checkboxes.length; i++) {
    checkboxes[i].checked = source.checked;
  }
}

function verify_input_data_before_post() {
    if (!verify_single_row("design_file")) {
        return "Please supply a design file or IUPAC and submit again"
    }
    if (!verify_single_row("config_file")) {
        return "Please supply a config file and submit again"

    }
        if (!(verify_single_row("fastq_file") || mounted_file_was_supplied('csv_file_server'))) {
        return "Please supply a fastq file and submit again"
     }

    return true;
}

function verify_single_row(row_name){
    var mapping = mapp_row_name_to_value[row_name]
    result = url_was_supplied(mapping["url"])|| upload_file_was_supplied(mapping["upload"])|| mounted_file_was_supplied(mapping["mounted"]);
    //if we are checking design file check if IUPAC was supplied
    if(row_name =="design_file"){
        return result || iupac_was_supplied();
    }
    return result;
}
//checked
function iupac_was_supplied() {
    var textbox_value = document.getElementById("iupac-textbox").value;
    var checkbox_value = document.getElementById("iupac-checkbox").checked;
        return (!["","IUPAC"].includes(textbox_value))&& checkbox_value;
}
//checked and works
function url_was_supplied(name){
    var textbox = document.getElementsByName(name)[0];
    console.log("inside url was supplied:"+ textbox);
    return (textbox.value != "")

}


//this func sends a get request to server and update reads files accordingly
function toggle_between_csv_fastq() {

    if($("input[type='checkbox'][name='after_matching']").is(':checked')) {
        uncheck_files("fastq_file_server")
        //uncheck choose all button
        document.getElementById('choose_all').checked = false;
        $("#csv_file_server").show()
        $("#fastq_file_server").hide()

    }else {
        uncheck_files("csv_file_server")
        $("#csv_file_server").hide()
        $("#fastq_file_server").show()
    }

}

function set_reads(extension){
    //change
    var response =  $.get('/extension/' + extension,
    (data)=> {
        var files = JSON.parse(data)
        update_reads(files)
    });
}
function uncheck_files(name) {

  checkboxes = document.getElementsByName(name);
  for(var i=0; i < checkboxes.length; i++) {
    checkboxes[i].checked = false;
  }

    // $('input[type='checkbox'=]')
    
}

//checked and works
function mounted_file_was_supplied(name){
    console.log("Hi from mounted file:" + name);
    checkboxes = document.getElementsByName(name);
    for(var i=0; i < checkboxes.length; i++) {
        if(checkboxes[i].checked && checkboxes[i].id !== "iupac-checkbox")
            return true;
    }
    return false;
}
//checked and works
function upload_file_was_supplied(name){
    return document.getElementsByName(name)[0].value!== ""
}


// verify that the user only chose one from checkbox except in fastq files

$('input:checkbox').on('change', function(){
    console.log("inside on change");
    if(treat_as_radio_button.includes(this.name)) {
        $('input[name="' + this.name + '"]').not(this).prop('checked', false);
        if(this.name == 'design_file_server'){
            $('#iupac-textbox').attr('disabled', true);

        }
    }
    //for IUPAC File
    if(this.id == 'iupac-checkbox'){
        if($('#iupac-checkbox').is(':checked')){
            console.log('Hi 1');
            $('#iupac-textbox').attr('disabled', false);
            }
         else{
             console.log('Hi 2');
             console.log($('#iupac-textbox'));
             $('#iupac-textbox').attr('disabled', true);
            }
        }

});


$(document).ready(function () {
    $('[data-toggle="tooltip"]').tooltip();
    //initially hide csv files
    $("#csv_file_server").hide();

    var design_file_list = document.getElementsByName('design_file_server');
    console.log('design file list:', design_file_list);
    if(design_file_list.length == 1){
         $('#iupac-checkbox').attr('checked', true);
    }
    else {

          for(var i=0; i < design_file_list.length; i++) {
              console.log("INSIDE FOR");
                if(design_file_list[i].id != 'iupac-checkbox'){
                    console.log('INSIDE IF');
                    design_file_list[i].checked = true;
                    break;
            }
  }
    }


    $('#data_submit').submit(function (event) {
        var token = new Date().getTime(); //use the current timestamp as the token value
        $('#time_stamp').val(token);
        var result = verify_input_data_before_post()
        if(result!==true) {
            alert(result)
            event.preventDefault();
            return;
        }

        $("#data_submit").hide();
        $("#uploading_text").text(" The files are being uploaded to the server, please wait.");
        $.blockUI({ message: $("#uploading_div")});

	});

    document.getElementById('choose_all').checked = true;
    toggle(document.getElementById('choose_all'));
  });
